import pandas as pd
import os
import EpiClockNBL.util as nbl_util
nbl_consts = nbl_util.consts

proj_dir_TARGET = os.path.join(nbl_consts['official_indir'], 'TARGET')
proj_dir_Henrich = os.path.join(nbl_consts['official_indir'], 'Henrich')

def pipeline():

    clinical = {}
    clinical['TARGET'] = pd.read_table(os.path.join(proj_dir_TARGET, 'clinical.processed.withAges.tsv'))
    clinical['Henrich'] = pd.read_table(os.path.join(proj_dir_Henrich, 'clinical.processed.withAges.tsv'))

    # Features to summarize in each cohort
    summary_cat_features = {
        'TARGET': [
            'classification_of_tumor', 'primary_diagnosis', 'morphology',
            'inss_stage', 'cog_neuroblastoma_risk_group', 'MYCN status', 'primary_site',
            'race', 'sex_at_birth', 'ethnicity', 'vital_status'
        ],
        'Henrich': [
            'age at diagnosis', 'inss stage', 'current risk category', 'mycn status',
            '1p status', '11q status', '17q status'
        ]
    }
    summary_num_features = {
        'TARGET': [
            'Age', 'year_of_diagnosis', 'days_to_last_follow_up'
        ],
        'Henrich': []
    }
    cat_orders_dict = {
        # TARGET
        'inss_stage': ['Stage 1', 'Stage 2B', 'Stage 3', 'Stage 4S', 'Stage 4'],
        'cog_neuroblastoma_risk_group': ['Low risk', 'Intermediate risk', 'High risk'],
        'ethnicity': ['not hispanic or latino', 'hispanic or latino', 'Unknown'],
        # Henrich
        'age at diagnosis': ['Age < 1.5 years', 'Age >= 1.5 years'],
        'inss stage': ['Stage 1', 'Stage 2', 'Stage 3', 'Stage 4S', 'Stage 4'],
        'current risk category': ['Low risk', 'Intermediate risk', 'High risk']
    }
    replace_features_dict = {
        # TARGET
        'inss_stage': 'INSS stage',
        'cog_neuroblastoma_risk_group': 'COG neuroblastoma risk group',
        'MYCN status': 'MYCN status',
        'inpc_grade': 'INPC grade',
        'Age': 'Age, years',
        # 'Percent Tumor Min':'Percent tumor'
        # Henrich
        'inss stage': 'INSS stage',
        'mycn status': 'MYCN status'
    }

    def capFirstLetter(x):
        return x[0].upper() + x[1:]
    def replaceFeatures(feat_name):
        return replace_features_dict.get(feat_name, feat_name.replace('_', ' ').capitalize())
    def featNameRow(feat_name):
        return [replaceFeatures(feat_name), '']
    def featCountRows(ser):
        counts_df = ser.value_counts(dropna=False).to_frame()
        col_name = counts_df.index.name

        # Order the categories
        cat_order = cat_orders_dict.get(col_name, counts_df.index)
        assert len(set(counts_df.index).difference(set(cat_order))) == 0
        # cat_order = sorted(cat_order, key=lambda x:cat_orders_dict.index(x))
        cat_order = list(filter(lambda x:x in counts_df.index, cat_order))
        counts_df = counts_df.loc[cat_order].copy()
        
        pct_float = (counts_df['count'] / counts_df['count'].sum() * 100)
        pct_str = pct_float.apply(formatValue)
        # Number and percentage together: N (XX%)
        counts_df['count'] = counts_df['count'].astype(str) + ' (' + pct_str + '%)'
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

    # for name, clinical_tbl in zip(['all', 'omit_ganglios', 'ganglios', 'lump_pure'], [clinical['TARGET'],
    #                                                                     clinical['TARGET'].loc[clinical['TARGET']['primary_diagnosis'] != 'Ganglioneuroblastoma'],
    #                                                                     clinical['TARGET'].loc[clinical['TARGET']['primary_diagnosis'] == 'Ganglioneuroblastoma'],
    #                                                                     clinical['TARGET'].loc[clinical['TARGET']['LUMP'] >= 0.7]
    #                                                                     ]):
    # Create a table for each cohort
    for cohort in ['TARGET', 'Henrich']:
        clinical_tbl = clinical[cohort]
        
        table_text_list = []
        table_text_list.append(
            ['Patients', f'{clinical_tbl.shape[0]} (100%)']
        )
        for num_feat in summary_num_features[cohort]:
            table_text_list.append(featNameRow(num_feat))
            table_text_list.extend(featDescribeRows(clinical_tbl[num_feat]))
        for cat_feat in summary_cat_features[cohort]:
            table_text_list.append(featNameRow(cat_feat))
            try:
                table_text_list.extend(featCountRows(clinical_tbl[cat_feat]))
            except:
                print(cat_feat)
                raise
        pd.DataFrame(data=table_text_list, columns=['Characteristic', 'Value'],
                    ).to_excel(f'{cohort}_patient_characteristics.xlsx', index=False)