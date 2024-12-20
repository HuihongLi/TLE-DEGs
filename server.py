import threading
import time
import requests
import os
import base64  
from flask import Flask
import pandas as pd
import numpy as np
import plotly.express as px
from dash import Dash, dcc, html, Input, Output, State
import dash_bootstrap_components as dbc

df1 = pd.read_csv("data/DEseq2.csv")
df2 = pd.read_csv("data/edgeR.csv")
df3 = pd.read_csv("data/limma.csv")
df4 = pd.read_csv("data/Wilcoxon.csv")

df1['cluster'] = 'DESeq2'
df2['cluster'] = 'edgeR'
df3['cluster'] = 'limma'
df4['cluster'] = 'Wilcoxon'

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H2("Threshold Selection", className="text-center mb-4"),
            dbc.Card([
                html.H4("DESeq2", className="card-title"),
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Minimum Absolute Fold Change (|log2FC|):"),
                        dbc.Input(id='fc_min_threshold_deseq2', type='number', value=1, step=0.5),
                    ], width=12),
                    dbc.Col([
                        dbc.Label("Maximum Absolute Fold Change (|log2FC|):"),
                        dbc.Input(id='fc_max_threshold_deseq2', type='number', value=5, step=0.5),
                    ], width=12),
                    dbc.Col([
                        dbc.Label("Adjusted P-value Threshold (padj):"),
                        dbc.Input(id='padj_threshold_deseq2', type='number', value=0.05, step=0.01),
                    ], width=12),
                ]),
            ], body=True, className='mb-3'),
            dbc.Card([
                html.H4("edgeR", className="card-title"),
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Minimum Absolute Fold Change (|log2FC|):"),
                        dbc.Input(id='fc_min_threshold_edger', type='number', value=1, step=0.5),
                    ], width=12),
                    dbc.Col([
                        dbc.Label("Maximum Absolute Fold Change (|log2FC|):"),
                        dbc.Input(id='fc_max_threshold_edger', type='number', value=5, step=0.5),
                    ], width=12),
                    dbc.Col([
                        dbc.Label("Adjusted P-value Threshold (padj):"),
                        dbc.Input(id='padj_threshold_edger', type='number', value=0.05, step=0.01),
                    ], width=12),
                ]),
            ], body=True, className='mb-3'),
            dbc.Card([
                html.H4("limma", className="card-title"),
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Minimum Absolute Fold Change (|log2FC|):"),
                        dbc.Input(id='fc_min_threshold_limma', type='number', value=1, step=0.5),
                    ], width=12),
                    dbc.Col([
                        dbc.Label("Maximum Absolute Fold Change (|log2FC|):"),
                        dbc.Input(id='fc_max_threshold_limma', type='number', value=5, step=0.5),
                    ], width=12),
                    dbc.Col([
                        dbc.Label("Adjusted P-value Threshold (padj):"),
                        dbc.Input(id='padj_threshold_limma', type='number', value=0.05, step=0.01),
                    ], width=12),
                ]),
            ], body=True, className='mb-3'),
            dbc.Card([
                html.H4("Wilcoxon", className="card-title"),
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Minimum Absolute Fold Change (|log2FC|):"),
                        dbc.Input(id='fc_min_threshold_wilcoxon', type='number', value=1, step=0.5),
                    ], width=12),
                    dbc.Col([
                        dbc.Label("Maximum Absolute Fold Change (|log2FC|):"),
                        dbc.Input(id='fc_max_threshold_wilcoxon', type='number', value=5, step=0.5),
                    ], width=12),
                    dbc.Col([
                        dbc.Label("Adjusted P-value Threshold (padj):"),
                        dbc.Input(id='padj_threshold_wilcoxon', type='number', value=0.05, step=0.01),
                    ], width=12),
                ]),
            ], body=True),
            dbc.Button(
                'Get List',
                id='download_button',
                n_clicks=0,
                color='primary',
                className='mt-4',
                style={'width': '100%'}
            ),
            dcc.Download(id='download-intersection'),
            dbc.Button(
                'Analysis',
                id='analysis_button',
                n_clicks=0,
                color='success',
                className='mt-2',
                style={'width': '100%'}
            ),
        ], width=2),
        dbc.Col([
            dcc.Graph(id='scatter_plot', style={'height': '85vh'}),
            html.Div(id='analysis_output', className='mt-4'),
        ], width=10)
    ], className="mt-4")
], fluid=True)

@app.callback(
    Output('scatter_plot', 'figure'),
    Input('fc_min_threshold_deseq2', 'value'),
    Input('fc_max_threshold_deseq2', 'value'),
    Input('padj_threshold_deseq2', 'value'),
    Input('fc_min_threshold_edger', 'value'),
    Input('fc_max_threshold_edger', 'value'),
    Input('padj_threshold_edger', 'value'),
    Input('fc_min_threshold_limma', 'value'),
    Input('fc_max_threshold_limma', 'value'),
    Input('padj_threshold_limma', 'value'),
    Input('fc_min_threshold_wilcoxon', 'value'),
    Input('fc_max_threshold_wilcoxon', 'value'),
    Input('padj_threshold_wilcoxon', 'value'),
)
def update_plot(fc_min_deseq2, fc_max_deseq2, padj_deseq2,
                fc_min_edger, fc_max_edger, padj_edger,
                fc_min_limma, fc_max_limma, padj_limma,
                fc_min_wilcoxon, fc_max_wilcoxon, padj_wilcoxon):
    df1_copy = df1.copy()
    abs_log2fc = df1_copy['log2FoldChange'].abs()
    df1_copy['Label'] = np.where(
        (abs_log2fc >= fc_min_deseq2) & (abs_log2fc <= fc_max_deseq2) & (df1_copy['padj'] < padj_deseq2),
        np.where(df1_copy['log2FoldChange'] > 0, 'Up', 'Down'),
        'Non'
    )
    df2_copy = df2.copy()
    abs_log2fc = df2_copy['log2FoldChange'].abs()
    df2_copy['Label'] = np.where(
        (abs_log2fc >= fc_min_edger) & (abs_log2fc <= fc_max_edger) & (df2_copy['padj'] < padj_edger),
        np.where(df2_copy['log2FoldChange'] > 0, 'Up', 'Down'),
        'Non'
    )
    df3_copy = df3.copy()
    abs_log2fc = df3_copy['log2FoldChange'].abs()
    df3_copy['Label'] = np.where(
        (abs_log2fc >= fc_min_limma) & (abs_log2fc <= fc_max_limma) & (df3_copy['padj'] < padj_limma),
        np.where(df3_copy['log2FoldChange'] > 0, 'Up', 'Down'),
        'Non'
    )
    df4_copy = df4.copy()
    abs_log2fc = df4_copy['log2FoldChange'].abs()
    df4_copy['Label'] = np.where(
        (abs_log2fc >= fc_min_wilcoxon) & (abs_log2fc <= fc_max_wilcoxon) & (df4_copy['padj'] < padj_wilcoxon),
        np.where(df4_copy['log2FoldChange'] > 0, 'Up', 'Down'),
        'Non'
    )

    df = pd.concat([df1_copy, df2_copy, df3_copy, df4_copy], axis=0)
    df['neg_log10_padj'] = -np.log10(df['padj'].replace(0, np.nan))
    min_size = 5
    max_size = 100
    df['-log10(padj)'] = (df['neg_log10_padj'] - df['neg_log10_padj'].min()) / (
        df['neg_log10_padj'].max() - df['neg_log10_padj'].min()
    ) * (max_size - min_size) + min_size
    df = df.dropna(subset=['neg_log10_padj'])
    cluster_order = ['DESeq2', 'edgeR', 'limma', 'Wilcoxon']
    df['x_jitter'] = df['cluster'].map({c: i for i, c in enumerate(cluster_order)}) + np.random.uniform(-0.2, 0.2, len(df))

    fig = px.scatter(
        df,
        x='x_jitter',
        y='log2FoldChange',
        color='Label',
        size='-log10(padj)',
        hover_data={
            'Gene': df['gene'],       
            'log2FoldChange': True,
            'padj': True,
            '-log10(padj)': True,  
            'x_jitter': False,       
            "Label": False
        },
        color_discrete_map={'Up': 'red', 'Down': 'blue', 'Non': 'grey'},
        title="Epilepsy Model Biomarker Identification",
        labels={"log2FoldChange": "log2FoldChange", '-log10(padj)': 'PointSize -log10(padj)'},
    )
    fig.update_traces(marker=dict(opacity=0.7))
    fig.update_layout(
        xaxis=dict(
            tickmode='array',
            tickvals=list(range(len(cluster_order))),
            ticktext=cluster_order,
            title='Method',
            showgrid=False  
        ),
        yaxis=dict(
            title='log2FoldChange',
            zeroline=True,
            zerolinewidth=1.5,
            zerolinecolor='black',
        ),
        shapes=[
            dict(
                type='line',
                x0=i,
                x1=i,
                y0=df['log2FoldChange'].min(),
                y1=df['log2FoldChange'].max(),
                line=dict(color="LightGrey", width=1)
            ) for i in range(len(cluster_order))
        ],
        template='plotly_white',
        legend_title_text='Gene Regulation'
    )

    return fig

@app.callback(
    Output('download-intersection', 'data'),
    Input('download_button', 'n_clicks'),
    State('fc_min_threshold_deseq2', 'value'),
    State('fc_max_threshold_deseq2', 'value'),
    State('padj_threshold_deseq2', 'value'),
    State('fc_min_threshold_edger', 'value'),
    State('fc_max_threshold_edger', 'value'),
    State('padj_threshold_edger', 'value'),
    State('fc_min_threshold_limma', 'value'),
    State('fc_max_threshold_limma', 'value'),
    State('padj_threshold_limma', 'value'),
    State('fc_min_threshold_wilcoxon', 'value'),
    State('fc_max_threshold_wilcoxon', 'value'),
    State('padj_threshold_wilcoxon', 'value'),
    prevent_initial_call=True
)

def download_intersection(n_clicks,
                          fc_min_deseq2, fc_max_deseq2, padj_deseq2,
                          fc_min_edger, fc_max_edger, padj_edger,
                          fc_min_limma, fc_max_limma, padj_limma,
                          fc_min_wilcoxon, fc_max_wilcoxon, padj_wilcoxon):
    intersect_genes = get_intersect_genes(fc_min_deseq2, fc_max_deseq2, padj_deseq2,
                                          fc_min_edger, fc_max_edger, padj_edger,
                                          fc_min_limma, fc_max_limma, padj_limma,
                                          fc_min_wilcoxon, fc_max_wilcoxon, padj_wilcoxon)

    if len(intersect_genes) > 0:
        genes_text = '\n'.join(sorted(intersect_genes))
        return dict(content=genes_text, filename="intersection_genes.txt")
    else:
        return dict(content="No genes found in intersection.", filename="intersection_genes.txt")

def get_intersect_genes(fc_min_deseq2, fc_max_deseq2, padj_deseq2,
                        fc_min_edger, fc_max_edger, padj_edger,
                        fc_min_limma, fc_max_limma, padj_limma,
                        fc_min_wilcoxon, fc_max_wilcoxon, padj_wilcoxon):
    df1_copy = df1.copy()
    abs_log2fc = df1_copy['log2FoldChange'].abs()
    df1_copy['Label'] = np.where(
        (abs_log2fc >= fc_min_deseq2) & (abs_log2fc <= fc_max_deseq2) & (df1_copy['padj'] < padj_deseq2),
        np.where(df1_copy['log2FoldChange'] > 0, 'Up', 'Down'),
        'Non'
    )
    degs_deseq2 = df1_copy[df1_copy['Label'] != 'Non']['gene'].unique()
    df2_copy = df2.copy()
    abs_log2fc = df2_copy['log2FoldChange'].abs()
    df2_copy['Label'] = np.where(
        (abs_log2fc >= fc_min_edger) & (abs_log2fc <= fc_max_edger) & (df2_copy['padj'] < padj_edger),
        np.where(df2_copy['log2FoldChange'] > 0, 'Up', 'Down'),
        'Non'
    )
    degs_edger = df2_copy[df2_copy['Label'] != 'Non']['gene'].unique()
    df3_copy = df3.copy()
    abs_log2fc = df3_copy['log2FoldChange'].abs()
    df3_copy['Label'] = np.where(
        (abs_log2fc >= fc_min_limma) & (abs_log2fc <= fc_max_limma) & (df3_copy['padj'] < padj_limma),
        np.where(df3_copy['log2FoldChange'] > 0, 'Up', 'Down'),
        'Non'
    )
    degs_limma = df3_copy[df3_copy['Label'] != 'Non']['gene'].unique()
    df4_copy = df4.copy()
    abs_log2fc = df4_copy['log2FoldChange'].abs()
    df4_copy['Label'] = np.where(
        (abs_log2fc >= fc_min_wilcoxon) & (abs_log2fc <= fc_max_wilcoxon) & (df4_copy['padj'] < padj_wilcoxon),
        np.where(df4_copy['log2FoldChange'] > 0, 'Up', 'Down'),
        'Non'
    )
    degs_wilcoxon = df4_copy[df4_copy['Label'] != 'Non']['gene'].unique()
    intersect_genes = set(degs_deseq2) & set(degs_edger) & set(degs_limma) & set(degs_wilcoxon)
    return intersect_genes

@app.callback(
    Output('analysis_output', 'children'),
    Input('analysis_button', 'n_clicks'),
    State('fc_min_threshold_deseq2', 'value'),
    State('fc_max_threshold_deseq2', 'value'),
    State('padj_threshold_deseq2', 'value'),
    State('fc_min_threshold_edger', 'value'),
    State('fc_max_threshold_edger', 'value'),
    State('padj_threshold_edger', 'value'),
    State('fc_min_threshold_limma', 'value'),
    State('fc_max_threshold_limma', 'value'),
    State('padj_threshold_limma', 'value'),
    State('fc_min_threshold_wilcoxon', 'value'),
    State('fc_max_threshold_wilcoxon', 'value'),
    State('padj_threshold_wilcoxon', 'value'),
    prevent_initial_call=True
)

def perform_analysis(n_clicks,
                     fc_min_deseq2, fc_max_deseq2, padj_deseq2,
                     fc_min_edger, fc_max_edger, padj_edger,
                     fc_min_limma, fc_max_limma, padj_limma,
                     fc_min_wilcoxon, fc_max_wilcoxon, padj_wilcoxon):
    intersect_genes = get_intersect_genes(fc_min_deseq2, fc_max_deseq2, padj_deseq2,
                                          fc_min_edger, fc_max_edger, padj_edger,
                                          fc_min_limma, fc_max_limma, padj_limma,
                                          fc_min_wilcoxon, fc_max_wilcoxon, padj_wilcoxon)

    if len(intersect_genes) == 0:
        return html.Div("No genes found in intersection.")
    genes_list = list(intersect_genes)
    species_id = 10090  
    string_api_url = "https://string-db.org/api"
    output_format = "image"
    method = "enrichmentfigure"
    request_url = "/".join([string_api_url, output_format, method])
    params = {
        "identifiers": "%0d".join(genes_list),
        "species": species_id,
        "caller_identity": "your_app_name",  # Replace with your app name or email
    }
    response = requests.post(request_url, data=params)
    if response.status_code != 200:
        return html.Div("Failed to retrieve enrichment figure from STRING API.")
    enrichment_image = "data:image/png;base64,{}".format(
        base64.b64encode(response.content).decode()
    )
    return html.Div([
        html.H4("Enrichment Analysis:"),
        html.Img(src=enrichment_image, style={'width': '100%', 'height': 'auto'}),
    ])

def keep_alive():
    """Ping the app periodically to prevent it from spinning down."""
    port = os.environ.get('PORT', 10000)
    url = f"http://localhost:{port}"
    while True:
        try:
            response = requests.get(url)
            if response.status_code == 200:
                print("Keep-alive request successful.")
            else:
                print(f"Keep-alive request failed with status {response.status_code}.")
        except Exception as e:
            print(f"Keep-alive request error: {e}")
        time.sleep(300)  
        
if __name__ == '__main__':
    threading.Thread(target=keep_alive, daemon=True).start()
    port = int(os.environ.get("PORT", 10000)) 
    app.run_server(debug=False, host='0.0.0.0', port=port)