"""Regenerate tutorial figures from committed simulation results (no simulations run).

Requires pandas and plotly. Run from any directory; --preview-dir additionally
exports PNGs for visual QA when Kaleido and Chrome are available.
"""
import argparse
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT = Path(__file__).resolve().parents[1] / 'docs' / 'parameter_selection'
TECHS = [('illumina_HS25', 'Illumina HiSeq2500'), ('ont_Q20', 'ONT Q20 model'), ('ont_Q25', 'ONT Q25 model')]
SETTINGS = [(0.05, 'p = 0.05', '#e67e22'), (1e-6, 'p = 1e-6', '#2980b9')]


def save(fig, name, height, preview, bottom=100):
    fig.update_layout(template='plotly_white', height=height,
                      font=dict(family='Arial, sans-serif', color='#2F3542', size=13),
                      legend=dict(orientation='h', y=-0.24, x=0.5, xanchor='center'),
                      margin=dict(t=65, b=bottom, l=60, r=20))
    fig.write_html(ROOT / f'{name}.html', include_plotlyjs='cdn', full_html=True,
                   div_id=name, config=dict(responsive=True, displaylogo=False))
    if preview:
        preview.mkdir(parents=True, exist_ok=True)
        fig.write_image(preview / f'{name}.png', width=800, height=height)


def main(preview=None):
    pairs = pd.read_csv(ROOT / 'data/simulation_pairs.csv')
    mixtures = pd.read_csv(ROOT / 'data/simulation_mixtures.csv')
    assert not pairs.duplicated(['technology', 'nominal_depth', 'error_rate', 'min_freq', 'p_threshold', 'kind', 'seed']).any()
    assert (pairs.recovered_differences + pairs.masked_differences + pairs.uncallable_differences == pairs.true_differences).all()
    pure = pairs[(pairs.kind == 'different') & (pairs.min_freq == .01)]
    fig = make_subplots(rows=1, cols=3, subplot_titles=[name for _, name in TECHS], shared_yaxes=True)
    for col, (tech, _) in enumerate(TECHS, 1):
        for p, label, color in SETTINGS:
            d = pure[(pure.technology == tech) & (pure.p_threshold == p)]
            g = d.groupby('nominal_depth').agg(mean=('masked_differences', 'mean'),
                    low=('masked_differences', 'min'), high=('masked_differences', 'max'),
                    n=('seed', 'size'), missing=('uncallable_differences', 'mean')).reset_index()
            fig.add_trace(go.Scatter(x=g.nominal_depth, y=g['mean'], name=label,
                legendgroup=label, showlegend=col == 1, mode='lines+markers', line=dict(color=color),
                error_y=dict(type='data', symmetric=False, array=g.high-g['mean'], arrayminus=g['mean']-g.low),
                customdata=g[['n', 'missing']],
                hovertemplate='Depth=%{x}x<br>Masked=%{y}/180<br>Replicates=%{customdata[0]}<br>Uncallable=%{customdata[1]}<extra>%{fullData.name}</extra>'), row=1, col=col)
    fig.update_xaxes(type='log', range=[0.9, 3.1], tickvals=[10, 30, 100, 300, 1000], title_text='Nominal depth (x)')
    fig.update_yaxes(range=[-1, 30])
    fig.update_yaxes(title_text='Masked true differences / 180', row=1, col=1)
    fig.update_layout(title=dict(text='Figure 4 · Recovered differences', x=.04, y=.99, yanchor='top', font=dict(size=16)))
    save(fig, 'fig4_simulation_masking', 455, preview)

    d = pairs[(pairs.kind == 'different') & (pairs.technology == 'ont_Q20') &
              (pairs.nominal_depth == 1000) & (pairs.p_threshold == .05)].sort_values('min_freq')
    assert d.masked_differences.tolist() == [94, 3]
    fig = go.Figure(go.Bar(x=['No frequency floor', '1% frequency floor'],
        y=d.masked_differences, marker_color=['#e67e22', '#2980b9'],
        text=[f'{n}/180 masked' for n in d.masked_differences], textposition='outside',
        customdata=d[['popani']], hovertemplate='%{x}<br>Masked=%{y}/180<br>popANI=%{customdata[0]:.4f}%<extra></extra>'))
    fig.update_yaxes(title_text='Masked true differences / 180', range=[0, 110])
    fig.update_layout(title=dict(text='Figure 5 · Frequency cutoff', x=.04, y=.99, yanchor='top', font=dict(size=16)), showlegend=False)
    save(fig, 'fig5_simulation_frequency_floor', 395, preview, bottom=75)

    fig = make_subplots(rows=1, cols=3, subplot_titles=[name for _, name in TECHS], shared_yaxes=True)
    for col, (tech, _) in enumerate(TECHS, 1):
        for p, label, color in SETTINGS:
            d = mixtures[(mixtures.technology == tech) & (mixtures.p_threshold == p)].sort_values('minor_fraction')
            fig.add_trace(go.Bar(x=[f'{v:.0%}' for v in d.minor_fraction],
                y=d.retained_minor_sites/180*100, name=label, legendgroup=label,
                showlegend=col == 1, marker_color=color, customdata=d[['retained_minor_sites', 'error_allele_sites']],
                hovertemplate='Minor strain=%{x}<br>Retained=%{customdata[0]}/180 (%{y:.1f}%)<br>Error-allele sites=%{customdata[1]}/17,820<extra>%{fullData.name}</extra>'), row=1, col=col)
    fig.update_yaxes(range=[0, 105])
    fig.update_yaxes(title_text='True minor-allele sites retained (%)', row=1, col=1)
    fig.update_xaxes(title_text='Minor-strain fraction')
    fig.update_layout(barmode='group', title=dict(text='Figure 6 · Rare-allele retention', x=.04, y=.99, yanchor='top', font=dict(size=16)))
    save(fig, 'fig6_simulation_minor_alleles', 455, preview)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preview-dir', type=Path)
    main(parser.parse_args().preview_dir)
