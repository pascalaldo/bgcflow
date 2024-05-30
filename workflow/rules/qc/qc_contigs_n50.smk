checkpoint qc:
    input:
        seqfu=expand("data/processed/{name}/tables/df_seqfu_stats.csv", name=PROJECT_IDS),
        checkm=expand("data/processed/{name}/tables/df_checkm_stats.csv", name=PROJECT_IDS),
    output: "data/processed/qc/qc_passed.csv",
    run:
        import pandas as pd
        df_seqfu = None
        for f in input.seqfu:
            df_s = pd.read_csv(f, sep=",", index_col=0, header=0)
            print(f)
            print(df_s)
            df_seqfu = pd.concat((df_seqfu, df_s))
        df_checkm = None
        for f in input.checkm:
            df_s = pd.read_csv(f, sep=",", index_col=0, header=0)
            print(f)
            print(df_s)
            df_checkm = pd.concat((df_checkm, df_s))
        print(df_seqfu)
        print(df_checkm)
        print("Final DF")
        df = df_seqfu.merge(df_checkm, how='outer', left_index=True, right_index=True)
        print(df)
        df["passed"] = ((df["N50"] > 50_000) & (df["# contigs"] < 200))
        df.to_csv(output, columns=["passed"], header=True, index=True, index_label="genome_id")
