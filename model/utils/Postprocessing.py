# Post processing for OSeMOSYS short

import pandas as pd
from .OSeMOSYS_PULP_functions import *
def postprocessing(res_df, df, sets_df, defaults_df):
    #Setup
    defaults_df = defaults_df
    roa = res_df[res_df['NAME']=='RateOfActivity'].copy()
    roa['YEAR'] = roa['YEAR'].astype(int)
    roa['MODE_OF_OPERATION'] = roa['MODE_OF_OPERATION'].astype(int)
    NSC = res_df[res_df['NAME']=='NewStorageCapacity'].copy()
    NSC['YEAR'] = NSC['YEAR'].astype(int)
    oar = df[df['PARAM']=='OutputActivityRatio']
    oar = oar[oar['VALUE'] >0].copy(deep=False)
    oar_f = oar[['FUEL', 'MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR', 'VALUE']].copy()
    oar_f['YEAR'] = oar_f['YEAR'].astype(int)
    oar_f['MODE_OF_OPERATION'] = oar_f['MODE_OF_OPERATION'].astype(int)
    #oar_f['TECHNOLOGY'] = oar_f['MODE_OF_OPERATION'].astype(str)
    iar = df[df['PARAM']=='InputActivityRatio']
    iar_f = oar[['FUEL', 'MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR', 'VALUE']].copy()
    iar_f['YEAR'] = iar_f['YEAR'].astype(int)
    iar_f['MODE_OF_OPERATION'] = iar_f['MODE_OF_OPERATION'].astype(int)
    ear = df[df['PARAM']=='EmissionActivityRatio']
    ear_f = ear[['EMISSION', 'MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR', 'VALUE']].copy()
    ear_f['YEAR'] = ear_f['YEAR'].astype(int)
    ear_f['MODE_OF_OPERATION'] = ear_f['MODE_OF_OPERATION'].astype(int)
    dr = df[df['PARAM']=='DiscountRateTech']
    dr_f = dr[['TECHNOLOGY', 'VALUE']].copy()
    
    #ep_f = df[df['PARAM']=='EmissionsPenalty'][['EMISSION', 'YEAR', 'VALUE']].copy()
    #if (ep_f['VALUE']!=0).any():
    #    ep_f['YEAR'] = ep_f['YEAR'].astype(int)
    #else:
    #    pass
    
    ys = df[df['PARAM']=='YearSplit']
    ys_f = ys[['TIMESLICE', 'YEAR', 'VALUE']].copy()
    ys_f['YEAR'] = ys_f['YEAR'].astype(int)
    omoo = df[df['PARAM']=='OutputModeofoperation']
    omoo_f = omoo[['MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR', 'VALUE']].copy()
    omoo_f['MODE_OF_OPERATION'] = omoo_f['MODE_OF_OPERATION'].astype(int)
    NSC = res_df[res_df['NAME']=='NewStorageCapacity'].copy()
    NSC['YEAR'] = NSC['YEAR'].astype(int)

    if len(sets_df['STORAGE']) != 0:
        CCS = df[df['PARAM']=='CapitalCostStorage']
        CCS = CCS[['REGION', 'STORAGE', 'YEAR', 'VALUE']].copy()
        drs = df[df['PARAM']=='DiscountRateSto']
        drs_f = drs[['STORAGE', 'VALUE']].copy()
    #Production By technology
    df_merge = pd.merge(roa, oar_f, on=['MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR'])
    df_merge = pd.merge(df_merge, ys_f, on=['TIMESLICE', 'YEAR'])
    namelist = []
    for i in range(0,len(df_merge)):
        namelist.append('ProductionByTechnology')
    df_PBT = df_merge
    df_PBT['VALUE_r'] = df_PBT['VALUE'] * df_PBT['VALUE_x'] * df_PBT['VALUE_y']
    df_PBT['NAME'] = namelist
    df_PBT = df_PBT.drop(['VALUE_x', 'VALUE', 'VALUE_y', 'FUEL_x'], axis=1)
    df_PBT.rename(columns = {'VALUE_r':'VALUE', 'FUEL_y': 'FUEL' }, inplace = True)
    res_df = pd.concat([res_df, df_PBT], ignore_index=True, sort=False)

    #ProductionbyTechnologyAnnual

    df2 = df_PBT.groupby(['YEAR', 'TECHNOLOGY', 'NAME']).sum()
    df2=df2.reset_index()
    namelist = []
    for i in range(0, len(df2)):
        namelist.append('ProductionByTechnologyAnnual')
    df2['NAME'] = namelist
    res_df = pd.concat([res_df, df2], ignore_index=True, sort=False)
    
    df_merge = pd.merge(roa, iar_f, on=['MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR'])
    df_merge = pd.merge(df_merge, ys_f, on=['TIMESLICE', 'YEAR'])

    #Use By technology
    namelist = []
    for i in range(0,len(df_merge)):
        namelist.append('UseByTechnology')
    df_UBT = df_merge
    df_UBT['VALUE_r'] = df_UBT['VALUE'] * df_UBT['VALUE_x'] * df_UBT['VALUE_y']
    df_UBT['NAME'] = namelist
    df_UBT = df_UBT.drop(['VALUE_x', 'VALUE', 'VALUE_y', 'FUEL_x'], axis=1)
    df_UBT.rename(columns = {'VALUE_r':'VALUE', 'FUEL_y': 'FUEL' }, inplace = True)
    res_df = pd.concat([res_df, df_UBT], ignore_index=True, sort=False)

    #UsebyTechnologyAnnual

    df2 = df_UBT.groupby(['YEAR', 'TECHNOLOGY', 'NAME']).sum()
    df2=df2.reset_index()
    namelist = []
    for i in range(0, len(df2)):
        namelist.append('UseByTechnologyAnnual')
    df2['NAME'] = namelist
    res_df = pd.concat([res_df, df2], ignore_index=True, sort=False)
    
    #AccumulatedNewCapacity
    newcap = res_df[res_df['NAME'] == 'NewCapacity']
    oplife = df[df['PARAM'] == 'OperationalLife'].copy()
    oplife = oplife[['PARAM', 'VALUE', 'TECHNOLOGY']]
    df_merge = pd.merge(oplife, newcap, on=['TECHNOLOGY'])

    ACC_df = pd.DataFrame()
    techlist = df_merge['TECHNOLOGY'].unique()
    for i in techlist:
        dfint = df_merge[df_merge['TECHNOLOGY'] == str(i)].copy()
        corvallist = []
        vallist = dfint['VALUE_y'].to_list()
        oplife = dfint['VALUE_x'].unique()[0]
        yearlist = dfint['YEAR'].to_list()
        for y in dfint['YEAR']:
            sum = 0
            for i in range(0, len(vallist)):
                if int(y) - int(yearlist[i]) >=0 and int(y) - int(yearlist[i]) < oplife:
                    sum = sum + vallist[i]
            corvallist.append(sum)
    newcap = res_df[res_df['NAME'] == 'NewCapacity']
    oplife = df[df['PARAM'] == 'OperationalLife'].copy()
    oplife = oplife[['PARAM', 'VALUE', 'TECHNOLOGY']]
    df_merge = pd.merge(oplife, newcap, on=['TECHNOLOGY'])

    ACC_df = pd.DataFrame()
    techlist = df_merge['TECHNOLOGY'].unique()
    for i in techlist:
        dfint = df_merge[df_merge['TECHNOLOGY'] == str(i)].copy()
        corvallist = []
        vallist = dfint['VALUE_y'].to_list()
        oplife = dfint['VALUE_x'].unique()[0]
        yearlist = dfint['YEAR'].to_list()
        for y in dfint['YEAR']:
            sum = 0
            for i in range(0, len(vallist)):
                if int(y) - int(yearlist[i]) >=0 and int(y) - int(yearlist[i]) < oplife:
                    sum = sum + vallist[i]
            corvallist.append(sum)
        dfint['VALUE'] = corvallist
        ACC_df = pd.concat([ACC_df, dfint], ignore_index=True, sort=False)
    namelist = []
    for i in range (0, len(ACC_df)):
        namelist.append('AccumulatedNewCapacity')
    ACC_df['NAME'] = namelist
    ACC_df = ACC_df.drop(['VALUE_x', 'PARAM', 'VALUE_y'], axis=1)
    res_df = pd.concat([res_df, ACC_df], ignore_index=True, sort=False)

    #AnnualTechnologyEmission
    df_merge = pd.merge(roa, ear_f, on=['MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR'])
    df_merge = pd.merge(df_merge, ys_f, on=['TIMESLICE', 'YEAR'])
    df_EBT = df_merge
    df_EBT['VALUE_r'] = df_EBT['VALUE'] * df_EBT['VALUE_x'] * df_EBT['VALUE_y']
    df_EBT = df_EBT.drop(['VALUE_x', 'VALUE', 'VALUE_y', 'EMISSION_x'], axis=1)
    df_EBT.rename(columns = {'VALUE_r':'VALUE', 'EMISSION_y': 'EMISSION' }, inplace = True)
    df2 = df_EBT.groupby(['YEAR', 'TECHNOLOGY', 'NAME', 'EMISSION']).sum()
    df2=df2.reset_index()
    namelist = []
    for i in range(0, len(df2)):
        namelist.append('AnnualTechnologyEmission')
    df2['NAME'] = namelist
    res_df = pd.concat([res_df, df2], ignore_index=True, sort=False)

    #DiscountedTechnologyEmissionsPenalty
    #if len(ep_f) != 0:
    #    ate = res_df[res_df['NAME']=='AnnualTechnologyEmission'].copy()
    #    df_merge = pd.merge(ate, ep_f, on=['EMISSION', 'YEAR'])
    #    df_merge = pd.merge(df_merge, dr_f, on=['TECHNOLOGY'])
    #    df_DEP = df_merge
    #    df_DEP['VALUE_r'] = df_DEP['VALUE_x'] * df_DEP['VALUE_y'] * ( 1 / (1 + df_DEP['VALUE'])**(df_DEP['YEAR'] - min(df_DEP['YEAR'])))
    #    df_DEP = df_DEP.drop(['VALUE_x', 'VALUE', 'VALUE_y'], axis=1)
    #    df_DEP.rename(columns = {'VALUE_r':'VALUE' }, inplace = True)
    #    namelist = []
    #    for i in range(0, len(df_DEP)):
    #        namelist.append('DiscountedTechnologyEmissionsPenalty')
    #    df_DEP['NAME'] = namelist
    #    res_df = pd.concat([res_df, df_DEP], ignore_index=True, sort=False)

    #TotalCapacityAnnual

    # Extract Residual Capacity and New Capacity data
    rescap = df[df['PARAM']=='ResidualCapacity'][['TECHNOLOGY', 'YEAR', 'VALUE' ]].copy()
    newcap = res_df[res_df['NAME'] == 'NewCapacity'][['TECHNOLOGY', 'YEAR', 'VALUE']].copy()
    newcap['YEAR'] = newcap['YEAR'].astype(int)
    #merge
    df_merge = pd.merge(rescap, newcap, on=['TECHNOLOGY', 'YEAR'])
    #Calculate totCapAnn
    df_merge['VALUE'] = df_merge['VALUE_x'] + df_merge['VALUE_y']
    df_merge['NAME'] = 'TotalCapacityAnnual'
    #Remove columns
    df_merge = df_merge.drop(['VALUE_x', 'VALUE_y',], axis=1)
    res_df = pd.concat([res_df, df_merge], ignore_index=True, sort=False)


    #DiscountedCapitalInvestment
    capcost = df[df['PARAM']=='CapitalCost'][['TECHNOLOGY', 'YEAR', 'VALUE']]
    newcap = res_df[res_df['NAME'] == 'NewCapacity'][['TECHNOLOGY', 'YEAR', 'VALUE']]
    newcap['YEAR']=newcap['YEAR'].astype(int)
    
    df_merge = pd.merge(capcost, newcap, on=['TECHNOLOGY', 'YEAR'])
    df_merge['YEAR']=df_merge['YEAR'].astype(int)
    
    dr_f = discount_factor(df=df, sets_df=sets_df, defaults_df=defaults_df)[0]
    df_merge = pd.merge(df_merge, dr_f, on=['TECHNOLOGY', 'YEAR'])
    
    df_merge['VALUE_r'] = df_merge['VALUE_x'] * df_merge['VALUE_y'] * df_merge['VALUE'] 
    df_merge = df_merge.drop(['VALUE_x', 'VALUE', 'VALUE_y'], axis=1)
    df_merge.rename(columns = {'VALUE_r':'VALUE' }, inplace = True)
    namelist = []
    for i in range(0, len(df_merge)):
        namelist.append('DiscountedCapitalInvestment')
    df_merge['NAME'] = namelist
    res_df = pd.concat([res_df, df_merge], ignore_index=True, sort=False)

    #AnnualFixedOperatingCost
    fixcost = df[df['PARAM']=='FixedCost']
    fixcost = fixcost[['TECHNOLOGY', 'VALUE']].copy()
    newcap = res_df[res_df['NAME'] == 'NewCapacity']
    df_merge = pd.merge(fixcost, newcap, on=['TECHNOLOGY'])
    df_merge['YEAR'] = df_merge['YEAR'].astype(int)
    df_merge['VALUE_r'] = df_merge['VALUE_x'] * df_merge['VALUE_y']
    df_merge = df_merge.drop(['VALUE_x', 'VALUE_y'], axis=1)
    df_merge.rename(columns = {'VALUE_r':'VALUE' }, inplace = True)
    namelist = []
    for i in range(0, len(df_merge)):
        namelist.append('AnnualFixedOperatingCost')
    df_merge['NAME'] = namelist
    res_df = pd.concat([res_df, df_merge], ignore_index=True, sort=False)


    #AnnualVariableOperatingCost
    varcost = df[df['PARAM']=='VariableCost']
    varcost = varcost[['MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR', 'VALUE']].copy()
    varcost['YEAR'] = varcost['YEAR'].astype(int)
    varcost['MODE_OF_OPERATION'] = varcost['MODE_OF_OPERATION'].astype(int)
    df_merge = pd.merge(roa, varcost, on=['MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR'])
    df_merge = pd.merge(df_merge, ys_f, on=['TIMESLICE', 'YEAR'])
    df_merge['VALUE_r'] = df_merge['VALUE'] * df_merge['VALUE_x'] * df_merge['VALUE_y']
    df_merge['VALUE_r'] = df_merge['VALUE_r'].astype(float)
    df_merge = df_merge.drop(['VALUE_x', 'VALUE', 'VALUE_y'], axis=1)
    df_merge.rename(columns = {'VALUE_r':'VALUE'}, inplace = True)
    df_merge = df_merge.groupby(['YEAR', 'TECHNOLOGY', 'NAME']).sum()
    df_merge=df_merge.reset_index()
    namelist = []
    for i in range(0, len(df_merge)):
        namelist.append('AnnualVariableOperatingCost')
    df_merge['NAME'] = namelist
    res_df = pd.concat([res_df, df_merge], ignore_index=True, sort=False)
    

    #DiscountedOperatingCost
    avc = res_df[res_df['NAME']=='AnnualVariableOperatingCost'].copy()
    afc = res_df[res_df['NAME']=='AnnualFixedOperatingCost'].copy()
    afc = afc[['TECHNOLOGY', 'YEAR', 'VALUE']].copy()
    df_merge = pd.merge(afc, avc, on=['TECHNOLOGY', 'YEAR'])
    df_merge = pd.merge(df_merge, dr_f, on=['TECHNOLOGY', 'YEAR'])
    df_merge['VALUE_r'] = (df_merge['VALUE_x'] + df_merge['VALUE_y']) * ( 1 / (1 + df_merge['VALUE'])**(df_merge['YEAR'] - min(df_merge['YEAR'])))
    df_merge = df_merge.drop(['VALUE_x', 'VALUE', 'VALUE_y'], axis=1)
    df_merge.rename(columns = {'VALUE_r':'VALUE' }, inplace = True)
    namelist = []
    for i in range(0, len(df_merge)):
        namelist.append('DiscountedOperatingCost')
    df_merge['NAME'] = namelist
    res_df = pd.concat([res_df, df_merge], ignore_index=True, sort=False)

    #TotalTechnologyAnnualActivity
   
    df_merge = pd.merge(roa, omoo_f, on=['MODE_OF_OPERATION', 'TECHNOLOGY', 'YEAR'])
    df_merge = pd.merge(df_merge, ys_f, on=['TIMESLICE', 'YEAR'])
    df_merge['VALUE_r'] = df_merge['VALUE'] * df_merge['VALUE_x'] * df_merge['VALUE_y']
    df_merge['VALUE_r']  = df_merge['VALUE_r'].astype(float)
    df_merge = df_merge.drop(['VALUE_x', 'VALUE', 'VALUE_y'], axis=1)
    df_merge.rename(columns = {'VALUE_r':'VALUE'}, inplace = True)
    df_merge = df_merge.groupby(['YEAR', 'TECHNOLOGY', 'NAME']).sum()
    df_merge=df_merge.reset_index()
    namelist = []
    for i in range(0, len(df_merge)):
        namelist.append('TotalTechnologyAnnualActivity')
    df_merge['NAME'] = namelist
    res_df = pd.concat([res_df, df_merge], ignore_index=True, sort=False)

    #TotalTechnologyModelPeriodActivity
    TTA = res_df[res_df['NAME']=='TotalTechnologyAnnualActivity'].copy()
    TTA['YEAR'] = TTA['YEAR'].astype(str)
    TTA['VALUE'] = TTA['VALUE'].astype(float)
    TMA = TTA.groupby(['TECHNOLOGY', 'NAME']).sum()
    TMA=TMA.reset_index()
    namelist = []
    for i in range(0, len(TMA)):
        namelist.append('TotalTechnologyModelPeriodActivity')
    TMA['NAME'] = namelist
    res_df = pd.concat([res_df, TMA], ignore_index=True, sort=False)

    #DiscountedCapitalInvestmentStorage
    if len(res_df[res_df['NAME']=='NewStorageCapacity'])!= 0:
        
        if len(sets_df['STORAGE']) != 0:
            stolist = CCS['STORAGE'].to_list()
            vallist = CCS['VALUE'].to_list()
            reglist = CCS['REGION'].to_list()
            yearlist = CCS['YEAR'].to_list()
            sets_df1 = sets_df[sets_df['STORAGE']!= 'nan']
            for i in sets_df1['STORAGE']:
                if i not in CCS['STORAGE'].unique():
                    for j in  sets_df['YEAR']:
                        if j!= 'nan':
                            stolist.append(i)
                            vallist.append(defaults_df[defaults_df['PARAM'] == 'CapitalCostStorage']['VALUE'].item())
                            reglist.append(sets_df['REGION'][0])
                            yearlist.append(j)

            CCS_new = pd.DataFrame()
            CCS_new['REGION'] = reglist
            CCS_new['VALUE'] = vallist
            CCS_new['STORAGE'] = stolist    
            CCS_new['YEAR'] = yearlist
            CCS_new['YEAR'] = CCS_new['YEAR'].astype(int)
            df_merge = pd.merge(CCS_new, NSC, on=['STORAGE', 'YEAR'])
            
            stolist = drs_f['STORAGE'].to_list()
            vallist = drs_f['VALUE'].to_list()
            sets_df1 = sets_df[sets_df['STORAGE']!= 'nan']
            for i in sets_df1['STORAGE']:
                    if i not in drs_f['STORAGE'].unique():
                        stolist.append(i)
                        vallist.append(defaults_df[defaults_df['PARAM'] == 'DiscountRateSto']['VALUE'].item())

            drs_fnew = pd.DataFrame()
            drs_fnew['VALUE'] = vallist
            drs_fnew['STORAGE'] = stolist
            df_merge = pd.merge(df_merge, drs_fnew, on=['STORAGE'])
            
            df_merge['YEAR'] = df_merge['YEAR'].astype(int)
            df_merge['VALUE_r'] = df_merge['VALUE_x'] * df_merge['VALUE_y'] * ( 1 / (1 + df_merge['VALUE'])**(df_merge['YEAR'] - min(df_merge['YEAR'])))
            df_merge = df_merge.drop(['VALUE_x', 'VALUE', 'VALUE_y'], axis=1)
            df_merge.rename(columns = {'VALUE_r':'VALUE' }, inplace = True)
            namelist = []
            for i in range(0, len(df_merge)):
                namelist.append('DiscountedCapitalInvestmentStorage')
            df_merge['NAME'] = namelist
            res_df = pd.concat([res_df, df_merge], ignore_index=True, sort=False)

        #ModelPeriodEmissions
        ate = res_df[res_df['NAME']=='AnnualTechnologyEmission'].copy()
        ate['YEAR'] = ate['YEAR'].astype(str)
        ate['VALUE'] = ate['VALUE'].astype(float)
        mpe = ate.groupby(['EMISSION', 'NAME']).sum()
        mpe=mpe.reset_index()
        namelist = []
        for i in range(0, len(mpe)):
            namelist.append('ModelPeriodEmissions')
        mpe['NAME'] = namelist
        res_df = pd.concat([res_df, mpe], ignore_index=True, sort=False)
    else:
        pass
    return res_df
    


