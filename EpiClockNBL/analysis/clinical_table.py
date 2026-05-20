import pandas as pd
import os
import EpiClockNBL.util as nbl_util
nbl_consts = nbl_util.consts

proj_dir = os.path.join(nbl_consts['official_indir'], 'TARGET')

def pipeline(verbose=True):
        
    clinical = {}
    clinical['TARGET'] = pd.read_table(os.path.join(proj_dir, 'clinical.TARGET.annotated.withMYCN.withPercentTumor.tsv'))

    summary_cat_features = [
        'classification_of_tumor', 'primary_diagnosis', 'morphology',
        'inss_stage', 'cog_neuroblastoma_risk_group', 'MYCN status',
        'race', 'gender', 'ethnicity', 'vital_status'
    ]
    summary_num_features = [
        'Age', 'Percent Tumor Min'
    ]
    cat_orders_dict = {
        'inss_stage': ['Stage 1', 'Stage 2B', 'Stage 3', 'Stage 4S', 'Stage 4']
    }
    replace_features_dict = {
        'inss_stage': 'INSS stage',
        'cog_neuroblastoma_risk_group': 'COG neuroblastoma risk group',
        'MYCN status': 'MYCN status',
        'inpc_grade': 'INPC grade',
        'Age': 'Age, years',
        'Percent Tumor Min':'Percent tumor'
    }

    def capFirstLetter(x):
        return x[0].upper() + x[1:]
    def replaceFeatures(feat_name):
        return replace_features_dict.get(feat_name, feat_name.replace('_', ' ').capitalize())
    def featNameRow(feat_name):
        return [replaceFeatures(feat_name), '', '']
    def featCountRows(ser):
        counts_df = ser.value_counts(dropna=False).to_frame()
        col_name = counts_df.columns[0]

        # Order the categories
        cat_order = cat_orders_dict.get(col_name, counts_df.index)
        assert len(set(counts_df.index).difference(set(cat_order))) == 0
        # cat_order = sorted(cat_order, key=lambda x:cat_orders_dict.index(x))
        cat_order = list(filter(lambda x:x in counts_df.index, cat_order))
        counts_df = counts_df.loc[cat_order].copy()
        
        pct_float = (counts_df[col_name] / counts_df[col_name].sum() * 100)
        pct_str = pct_float.apply(formatValue)
        counts_df['pct'] = pct_str
        counts_df[col_name] = counts_df[col_name].astype(str)
        return counts_df.rename(index=lambda x:capFirstLetter(str(x))).rename(index=lambda x:'Missing' if x.lower()=='nan' else x).rename(index=lambda x:' '*20 + x).reset_index().values.tolist()
    def formatValue(x):
        if x < 1:
            return f'{x:.2g}'
        return format(x, '.1f')
    def featDescribeRows(ser):
        describe_df = ser.describe().apply(lambda x:' '*20 + formatValue(x))
        col_name = describe_df.name
        describe_df = describe_df.loc[['min', '25%', '50%', '75%', 'max']].rename({'min':'Min', 'max':'Max'}).copy()
        return describe_df.rename(index=lambda x:' '*20 + x).reset_index().values.tolist()

    for name, clinical_tbl in zip(['all', 'omit_ganglios', 'ganglios', 'lump_pure'], [clinical['TARGET'],
                                                                        clinical['TARGET'].loc[clinical['TARGET']['primary_diagnosis'] != 'Ganglioneuroblastoma'],
                                                                        clinical['TARGET'].loc[clinical['TARGET']['primary_diagnosis'] == 'Ganglioneuroblastoma'],
                                                                        clinical['TARGET'].loc[clinical['TARGET']['LUMP'] >= 0.7]
                                                                        ]):
        table_text_list = []
        table_text_list.append(
            ['Patients', clinical_tbl.shape[0], '100']
        )
        for num_feat in summary_num_features:
            table_text_list.append(featNameRow(num_feat))
            table_text_list.extend(featDescribeRows(clinical_tbl[num_feat]))
        for cat_feat in summary_cat_features:
            table_text_list.append(featNameRow(cat_feat))
            try:
                table_text_list.extend(featCountRows(clinical_tbl[cat_feat]))
            except:
                print(name)
                print(cat_feat)
                raise
        pd.DataFrame(data=table_text_list, columns=['Characteristic', 'Number', '%'],
                    ).to_excel(f'{name}_patient_characteristics.xlsx', index=False)