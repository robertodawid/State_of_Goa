def rate_us_ts(df_ru):
    df_ru.dropna(axis=1)
    ts_list = df_ru['TIMESLICE'].tolist()
    ys_list = []

    for i in ts_list:
        if i == 1:
            ys_list.append(0.0685)
        elif i == 2:
            ys_list.append(0.1233)
        elif i == 3:
            ys_list.append(0.0548)
        elif i == 4:
            ys_list.append(0.0822)
        elif i == 5:
            ys_list.append(0.0696)
        elif i == 6:
            ys_list.append(0.1253)
        elif i == 7:
            ys_list.append(0.0557)
        elif i == 8:
            ys_list.append(0.0836)
        elif i == 9:
            ys_list.append(0.0702)
        elif i == 10:
            ys_list.append(0.1264)
        elif i == 11:
            ys_list.append(0.0562)
        elif i == 12:
            ys_list.append(0.0842)
            
    df_ru['Yearsplit'] = ys_list
    df_ru['Usebytech'] = df_ru['VALUE'] * df_ru['Yearsplit']
    df_ru.drop("Yearsplit", axis=1, inplace=True)
    df_ru.drop("VALUE", axis=1, inplace=True)
    df_ru.rename(columns={"Usebytech": "VALUE"}, inplace=True)

    df_ru = df_ru[
        [
        "NAME",
        "VALUE",
        "SCENARIO",
        "REGION",
        "REGION2",
        "DAYTYPE",
        "EMISSION",
        "FUEL",
        "DAILYTIMEBRACKET",
        "SEASON",
        "TIMESLICE",
        "MODE_OF_OPERATION",
        "STORAGE",
        "TECHNOLOGY",
        "YEAR",
    ]]
    return(df_ru)
