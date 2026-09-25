use super::{allowed, AnyResult, Args};
use rust_htslib::bam::{self, pileup, HeaderView, Read};
use rust_htslib::errors::{Error, Result};
use rust_htslib::tpool::ThreadPool;
use rust_htslib::{faidx, htslib};
use std::collections::HashMap;
use std::ffi::{c_void, CString};
use std::marker::PhantomData;
use std::slice;
use std::sync::{Arc, Mutex, OnceLock, Weak};

static REFERENCES: OnceLock<Mutex<HashMap<String, Weak<CString>>>> = OnceLock::new();

pub(super) struct BaqReader {
    inner: bam::IndexedReader,
    fasta: Option<faidx::Reader>,
    chromosome: String,
    reference: Option<Arc<CString>>,
    record: bam::Record,
    args: Args,
}

impl BaqReader {
    pub(super) fn new(args: &Args) -> Result<Self> {
        Ok(Self {
            inner: bam::IndexedReader::from_path(&args.bam)?,
            fasta: args
                .reference
                .as_deref()
                .map(faidx::Reader::from_path)
                .transpose()?,
            chromosome: String::new(),
            reference: None,
            record: bam::Record::new(),
            args: args.clone(),
        })
    }

    pub(super) fn fetch(&mut self, chrom: &str, start: i64, end: i64) -> Result<()> {
        if self.chromosome != chrom && self.fasta.is_some() {
            let fasta = self.fasta.as_ref().unwrap();
            let cache = REFERENCES.get_or_init(|| Mutex::new(HashMap::new()));
            let mut cache = cache.lock().map_err(|_| Error::FaidxBadSeqName)?;
            let sequence = if let Some(existing) = cache.get(chrom).and_then(Weak::upgrade) {
                existing
            } else {
                let len = fasta.fetch_seq_len(chrom) as usize;
                if len == 0 || len > i64::MAX as usize {
                    return Err(Error::FaidxBadSeqName);
                }
                let bases = fasta.fetch_seq(chrom, 0, len - 1)?;
                let value = Arc::new(CString::new(bases).map_err(|_| Error::FaidxBadSeqName)?);
                cache.insert(chrom.to_owned(), Arc::downgrade(&value));
                value
            };
            self.reference = Some(sequence);
        }
        self.chromosome = chrom.to_owned();
        self.inner.fetch((chrom, start, end))
    }

    pub(super) fn pileup_with_overlaps(&mut self, depth: u32) -> AnyResult<Mplp<'_>> {
        let mut data = [self as *mut Self as *mut c_void];
        let iter = unsafe { htslib::bam_mplp_init(1, Some(Self::next_read), data.as_mut_ptr()) };
        if iter.is_null() {
            return Err("could not initialize HTSlib mpileup".into());
        }
        if unsafe { htslib::bam_mplp_init_overlaps(iter) } < 0 {
            unsafe { htslib::bam_mplp_destroy(iter) };
            return Err("could not enable paired-read overlap handling".into());
        }
        unsafe { htslib::bam_mplp_set_maxcnt(iter, depth as i32) };
        Ok(Mplp {
            iter,
            _reader: PhantomData,
        })
    }

    extern "C" fn next_read(data: *mut c_void, dest: *mut htslib::bam1_t) -> i32 {
        let this = unsafe { &mut *(data as *mut Self) };
        loop {
            match this.inner.read(&mut this.record) {
                None => return -1,
                Some(Err(_)) => return -2,
                Some(Ok(())) => {
                    if this.record.flags() & 0x704 != 0 || !allowed(&this.record, &this.args) {
                        continue;
                    }
                    // Apply the same HTSlib BAQ adjustment used by mpileup when a FASTA is supplied.
                    unsafe {
                        if let Some(reference) = &this.reference {
                            let result = htslib::sam_prob_realn(
                                this.record.inner_mut(),
                                reference.as_ptr(),
                                reference.as_bytes().len() as i64,
                                htslib::htsRealnFlags_BAQ_APPLY as i32,
                            );
                            if result == -4 {
                                return -2;
                            }
                        }
                        return if htslib::bam_copy1(dest, this.record.inner()).is_null() {
                            -2
                        } else {
                            0
                        };
                    }
                }
            }
        }
    }
}

pub(super) struct Mplp<'a> {
    iter: htslib::bam_mplp_t,
    _reader: PhantomData<&'a mut BaqReader>,
}

pub(super) struct Position<'a> {
    tid: u32,
    pos: u32,
    alignments: &'a [htslib::bam_pileup1_t],
}

impl Position<'_> {
    pub(super) fn tid(&self) -> u32 {
        self.tid
    }
    pub(super) fn pos(&self) -> u32 {
        self.pos
    }
    pub(super) fn alignments(&self) -> impl Iterator<Item = bam::pileup::Alignment<'_>> {
        self.alignments.iter().map(bam::pileup::Alignment::new)
    }
}

impl Mplp<'_> {
    pub(super) fn next_position(&mut self) -> AnyResult<Option<Position<'_>>> {
        let mut tid = 0_i32;
        let mut pos = 0_i32;
        let mut depth = [0_i32; 1];
        let mut ptrs = [std::ptr::null(); 1];
        let result = unsafe {
            htslib::bam_mplp_auto(
                self.iter,
                &mut tid,
                &mut pos,
                depth.as_mut_ptr(),
                ptrs.as_mut_ptr(),
            )
        };
        if result < 0 {
            return Err("HTSlib mpileup failed".into());
        }
        if result == 0 {
            return Ok(None);
        }
        let alignments = unsafe { slice::from_raw_parts(ptrs[0], depth[0] as usize) };
        Ok(Some(Position {
            tid: tid as u32,
            pos: pos as u32,
            alignments,
        }))
    }
}

impl Drop for Mplp<'_> {
    fn drop(&mut self) {
        unsafe { htslib::bam_mplp_destroy(self.iter) }
    }
}

impl Read for BaqReader {
    fn read(&mut self, record: &mut bam::Record) -> Option<Result<()>> {
        loop {
            let next = self.inner.read(record);
            if !matches!(next, Some(Ok(()))) {
                return next;
            }
            if record.flags() & 0x704 != 0 || !allowed(record, &self.args) {
                continue;
            }
            if let Some(reference) = &self.reference {
                unsafe {
                    htslib::sam_prob_realn(
                        record.inner_mut(),
                        reference.as_ptr(),
                        reference.as_bytes().len() as i64,
                        htslib::htsRealnFlags_BAQ_APPLY as i32,
                    );
                }
            }
            return next;
        }
    }

    fn records(&mut self) -> bam::Records<'_, Self> {
        unreachable!("BaqReader is only used for pileup")
    }
    fn rc_records(&mut self) -> bam::RcRecords<'_, Self> {
        unreachable!("BaqReader is only used for pileup")
    }
    fn pileup(&mut self) -> pileup::Pileups<'_, Self> {
        let iter = unsafe {
            htslib::bam_plp_init(Some(Self::next_read), self as *mut Self as *mut c_void)
        };
        pileup::Pileups::new(self, iter)
    }
    fn htsfile(&self) -> *mut htslib::htsFile {
        self.inner.htsfile()
    }
    fn header(&self) -> &HeaderView {
        self.inner.header()
    }
    fn set_thread_pool(&mut self, pool: &ThreadPool) -> Result<()> {
        self.inner.set_thread_pool(pool)
    }
}
