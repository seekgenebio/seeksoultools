import os
import base64
import pandas as pd
import lzstring
import scipy.io as io
import plotly.express as px
import plotly.graph_objects as go

template_dir = os.path.join(os.path.dirname(__file__), 'template')
data = {}

def line_color(p):
    _color = ['#dbdbf4',
             '#c5c5f5',
             '#afaff6',
             '#9999f7',
             '#8383f8',
             '#6d6df9',
             '#5858fb',
             '#4141fc',
             '#2b2bfd',
             '#1616fe']
    if p==1:
        return '#0000ff'
    elif p==0:
        return '#dddddd'
    else:
        return _color[int(p/0.1)]



def diff_table(f, N=20):
    df = (pd.read_table(f)
         .assign(cluster=lambda df: df['cluster'].map(lambda x: f'cluster{x}'))
    )
    # table_id, classes
    fold_change_title = 'avg_logFC'
    if 'avg_log2FC' in df.columns:
        fold_change_title = 'avg_log2FC'
    df2 = (df.drop(columns=['p_val',  'pct.1',  'pct.2'])
        .groupby('cluster')
        .apply(lambda x: x.sort_values(['p_val_adj', fold_change_title], ascending=[True, False]).head(N))
        .reset_index(drop=True)
        .pivot_table(values=[fold_change_title, 'p_val_adj'], index=['Ensembl', 'gene'], columns='cluster')
        .swaplevel(axis=1)
        .sort_index(1, 0, key=lambda x: x.str.replace('cluster', '').astype(int))
    )
    # index_names
    return (df2.to_html(table_id='marker_table', classes='display', na_rep='-', index_names=False)
        .replace('border="1"', ''))


###zdef barcode_rank_data(countsFile, barcodesFile):
###z    df = pd.read_csv(countsFile, sep='\t')
###z    barcodes = pd.read_csv(barcodesFile, header=None, sep='\t')
###z    UMIcounts = df.groupby(['cellID'])['UMINum'].agg(UMIcounts=sum)
###z    UMIcounts = UMIcounts.sort_values(by='UMIcounts', ascending=False)
###z    UMIcounts = UMIcounts.reset_index()
    #UMIcounts['rank'] = UMIcounts.index + 1

    #max_cellRank = int(UMIcounts.loc[UMIcounts['cellID'].isin(barcodes[0]), 'rank'].iloc[-1])
###z    max_idx_cell = int(UMIcounts[UMIcounts['cellID'].isin(barcodes[0])].index[-1])
###z    data['cells_index'] = [0, max_idx_cell]  # rows=5, index=[0,1,2,3,4], 3 cell, [0,2] => [0, 1, 2] is cell
###z    data['cells_data'] = UMIcounts.UMIcounts[0: max_idx_cell + 1].to_list() # [0: 2+1] => [0, 1, 2]

    # start from the last cell to prevent line gap
###z    data['background_index'] = [max_idx_cell, int(UMIcounts.index[-1])]  # [2, 4] => [2, 3, 4] is background, start from 2 (cell) to prevent line gap
###z    data['background_data'] = UMIcounts.UMIcounts[max_idx_cell:].to_list() # [2: ] => [2, 3, 4]

def barcode_rank_plot_rna(cells_gz, barcodes_gz, mtx_gz):
    cells_df = pd.read_table(cells_gz, header=None, names=['barcode'])
    all_barcode_df = pd.read_table(barcodes_gz, header=None, names=['barcode'])
    cells_idx_df = all_barcode_df.loc[all_barcode_df.barcode.isin(cells_df.barcode),:]

    mat = io.mmread(mtx_gz)
    df = pd.DataFrame({"feature": mat.row, "barcode":mat.col, "umi": mat.data})
    umi_df = df.groupby(by="barcode").agg(umis=("umi", sum)).sort_values(by='umis', ascending=False)
    umi_df = umi_df.assign(is_cell = umi_df.index.isin(cells_idx_df.index)).reset_index()
    umi_df = umi_df.groupby(by=['umis', 'is_cell']).size().to_frame(name="count").reset_index().sort_values(by=['umis', 'is_cell'], ascending=False).reset_index(drop=True)
    return barcode_rank_plot(umi_df)



def barcode_rank_plot(umi_df):
    # continue
    idx0 = umi_df.loc[~umi_df['is_cell'],].index[0]
    idx1 = umi_df[::-1].loc[umi_df['is_cell'],].index[0]
    umi_df_head = umi_df.iloc[:idx0,]
    umi_df_tail = umi_df.iloc[idx1+1:,]
    umi_df_mix = umi_df.iloc[idx0:idx1+1,]

    # head, 1 trace
    res_list = []
    tmp = []
    n = 0
    if umi_df_head.shape[0]>0:
        for row in umi_df_head.itertuples(index=False):
        # for row in demo_df.itertuples(index=False):
            if (row.count==1):
                tmp.append((n, row.umis))
                n += row.count
            else:
                tmp.append((n, row.umis))
                tmp.append((n+row.count-1, row.umis))
                n += row.count
        res_list.append([tmp, f"100% Cells<br>{n}/{n}", 1])

    # mix
    tmp = []
    # row_num = umi_df_mix.shape[0]
    # break_num = int(np.log10(row_num)/20)
    break_num = n + 10
    step = 100
    counter = 0
    is_cell_num = 0
    barcode_num = 0
    for row in umi_df_mix.itertuples(index=False):
        if (row.count==1):
            tmp.append((n, row.umis))
            if row.is_cell:
                is_cell_num += row.count
            n += row.count
            barcode_num += row.count
        else:
            tmp.append((n, row.umis))
            tmp.append((n+row.count-1, row.umis))
            if row.is_cell:
                is_cell_num += row.count
            n += row.count
            barcode_num += row.count
        # counter += 1
        if  n >= break_num:
            # counter = 0
            res_list.append([tmp, f"{is_cell_num/barcode_num:.0%} Cells<br>{is_cell_num}/{barcode_num}", is_cell_num/barcode_num])
            tmp = []
            break_num = n + step
            step = step * 10
            is_cell_num = 0
            barcode_num = 0
            # print(break_num)
    else:
        if barcode_num !=0 :
            res_list.append([tmp, f"{is_cell_num/barcode_num:.0%} Cells<br>{is_cell_num}/{barcode_num}",  is_cell_num/barcode_num])

    # tail, 1 trace
    tmp = []
    for row in umi_df_tail.itertuples(index=False):
    # for row in demo_df.itertuples(index=False):
        if (row.count==1):
            tmp.append((n, row.umis, 1))
            n += row.count
        else:
            tmp.append((n, row.umis, 1))
            tmp.append((n+row.count-1, row.umis, 1))
            n += row.count
    res_list.append([tmp, "Background", 0])
    plot_data = []
    #print(res_list)
    for idx in range(len(res_list)):
        if idx > 0:
            plot_data.append({
                "x": [plot_data[-1]['x'][-1], ]  + [_[0] for _ in res_list[idx][0]],
                "y": [plot_data[-1]['y'][-1], ]  + [_[1] for _ in res_list[idx][0]],
                "text": res_list[idx][1],
                "hoverinfo": 'text',
                "line": dict(color=line_color(res_list[idx][-1])),
                "mode": 'lines',
                "legendgroup": 'Cells',
                "showlegend": False,
            })
        else:
            plot_data.append({
                "x": [_[0] for _ in res_list[idx][0]],
                "y": [_[1] for _ in res_list[idx][0]],
                "text": res_list[idx][1],
                "hoverinfo": 'text',
                "line": dict(color=line_color(res_list[idx][-1])),
                "mode": 'lines',
                "legendgroup": 'Cells',
                "showlegend": False,
            })

    plot_data[0]['name'] = 'Cells'
    plot_data[0]['showlegend'] = True
    plot_data[-1]['name'] = 'Background'
    plot_data[-1]['showlegend'] = True
    plot_data[-1]['legendgroup'] = 'Background'

    config = {
        "displayModeBar": True,
        "displaylogo": False,
        'modeBarButtonsToRemove': ['lasso','zoom','pan', 'zoomin', 'zoomout', 'autoscale', 'select2d'],
        'toImageButtonOptions': {
            'format': 'svg',
            'filename': 'custom_image',
            'height': 500,
            'width': 700,
            'scale': 1
          },
    }

    fig = go.Figure()
    for _ in plot_data:
        fig.add_trace(go.Scatter(**_))
    fig.update_layout(xaxis_type = "log", yaxis_type = "log")
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='#dee2e6')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='#dee2e6')
    fig.layout.height = 500
    fig.layout.width = 500
    fig.update_layout(
        #title='Barcode Rank',
        #title_x=0.5,
        title=None,
        xaxis=dict(title="Barcodes"),
        yaxis=dict(title="UMI counts"),
        xaxis_type = "log",
        yaxis_type = "log",
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        modebar_bgcolor='rgba(0,0,0,0)',
        modebar_color='#dee2e6'
    )
    # fig.show(config=config)
    json_dict = fig.to_plotly_json()
    json_dict['config'] = config
    return json_dict

def reduction_data(reduction_umi):
    df = pd.read_table(reduction_umi, index_col=0)
    df = df.drop(['orig.ident', 'nFeature_RNA', 'percent.mito'], axis=1)
    data['reduction'] = df.to_dict(orient='list')
    data['params'] = [ _ for _ in df.columns if _.startswith(('RNA_snn_res','seurat_clusters')) ]
    data['labs'] = df.columns[0:2].to_list()

def png2base64(f):
    with open(f, 'rb') as fh:
        return base64.b64encode(fh.read()).decode()

def report(samplename, outdir, **kwargs):
    import json
    from jinja2 import Environment, FileSystemLoader

    summary_file = os.path.join(outdir, f'{samplename}_summary.json')
    assert os.path.exists(summary_file), f'{summary_file} not found!'
    with open(summary_file) as fh:
        summary = json.load(fh)

    sequencing_table = {}
    sequencing_table['Number of Reads'] = f'{summary["stat"]["total"]:,}'
    sequencing_table['Valid Barcodes'] = f'{summary["stat"]["valid"]/summary["stat"]["total"]:.2%}'
    sequencing_table['Sequencing Saturation'] = f'{summary["cells"]["Sequencing Saturation"]:.2%}'
    del summary["cells"]["Sequencing Saturation"]

    if 'no_anchor' in summary["stat"]:
        sequencing_table['Without Anchor'] = f'{summary["stat"]["no_anchor"]:,}'

    if 'trimmed' in summary["stat"]:
        sequencing_table['Trimmed'] =  f'{summary["stat"]["trimmed"]:,}'
    if 'too_short' in summary["stat"]:
        sequencing_table['Too Short'] =  f'{summary["stat"]["too_short"]:,}'
    b_total_base = sum([sum(v) for v in summary["barcode_q"].values()])
    b30_base = sum([sum(v[30:]) for v in summary["barcode_q"].values()])
    sequencing_table['Q30 Bases in Barcode'] = f'{b30_base/b_total_base:.2%}'
    u_total_base = sum([sum(v) for v in summary["umi_q"].values()])
    u30_base = sum([sum(v[30:]) for v in summary["umi_q"].values()])
    sequencing_table['Q30 Bases in UMI'] = f'{u30_base/u_total_base:.2%}'

    mapping_table = {k: f'{v:.2%}' for k, v in summary["mapping"].items()}

    cells_table = dict([(k, f'{v:,}') if isinstance(v, int) else (k,f'{v:.2%}') for k,v in summary["cells"].items()])
    
    sample_table = {
        'Name': samplename, 
        'Description': '',
        'Transcriptome': summary["reference"],
        'Chemistry': summary["stat"].get("chemistry", "custom"),
        'Include introns': summary["include_introns"],
        'Seeksoul tools Version': summary["__version__"]
    }
    
    reduction_xls = os.path.join(outdir, 'step4', 'tsne_umi.xls')
    assert os.path.exists(reduction_xls), f'{reduction_xls} not found!'
    reduction_data(reduction_xls)

    count_xls = os.path.join(outdir, 'step3', 'counts.xls')
    barcodes_tsv = os.path.join(outdir, 'step3', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz')
    #barcode_rank_data(count_xls, barcodes_tsv)
    cells_gz = os.path.join(outdir, 'step3/filtered_feature_bc_matrix/barcodes.tsv.gz')
    barcodes_gz = os.path.join(outdir, 'step3/raw_feature_bc_matrix/barcodes.tsv.gz')
    mtx_gz = os.path.join(outdir, 'step3/raw_feature_bc_matrix/matrix.mtx.gz')
    barcod_rank_data = barcode_rank_plot_rna(cells_gz, barcodes_gz, mtx_gz)
    data["barcode_rank_data"]=barcod_rank_data

    f =  os.path.join(outdir, 'step4', 'FindAllMarkers.xls')
    marker_table = diff_table(f)

    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('base.html')
    with open(os.path.join(outdir, f'{samplename}_report.html'), 'w') as fh:
        rawdata = lzstring.LZString().compressToBase64(json.dumps(data))
        fh.write(template.render(
            logobase64 = png2base64(os.path.join(template_dir,'logo.png')),
            sequencing_table = sequencing_table,
            mapping_table = mapping_table,
            cells_table = cells_table,
            sample_table = sample_table,
            marker_table = marker_table,
            downsample = summary["downsample"],
            rawdata=rawdata)
        )

    header=('Samplename,Estimated_Number_of_Cells,Mean_Reads_per_Cell,Median_Genes_per_Cell,Number_of_Reads,'
            'Valid_Barcodes,Sequencing_Saturation,Reads_Mapped_Confidently_to_Genome,Fraction_Reads_in_Cells,'
            'Total_Genes_Detected,Median_UMI_Counts_per_Cell')

    summary_data = [
             samplename,
             cells_table['Estimated Number of Cells'],
#             sequencing_table['Number of Reads'],
             cells_table['Mean Reads per Cell'],
             cells_table['Median Genes per Cell'],
             sequencing_table['Number of Reads'],
             sequencing_table['Valid Barcodes'],
             sequencing_table['Sequencing Saturation'],
             mapping_table['Reads Mapped Confidently to Genome'],
             cells_table['Fraction Reads in Cells'],
             cells_table['Total Genes Detected'],
             cells_table['Median UMI Counts per Cell']
           ]

    with open(os.path.join(outdir, f'{samplename}_summary.csv'), 'w') as fh:
        fh.write(header + '\n')
        fh.write(','.join(str(_).replace(',', '') for _ in summary_data)+ '\n')

