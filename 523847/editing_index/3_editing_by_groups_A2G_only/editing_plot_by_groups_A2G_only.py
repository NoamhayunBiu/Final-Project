import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import mannwhitneyu, shapiro
from statsmodels.stats.multitest import multipletests

# -----------------------
# Helpers
# -----------------------
def rename_edit_col(df):
    # Unify the editing column name to EditingValue if needed
    if 'EditingValue' not in df.columns and 'A2GEditingIndex' in df.columns:
        df = df.rename(columns={'A2GEditingIndex':'EditingValue'})
    return df

def shapiro_by_group(df, group_cols=('Phase','treatment')):
    rows=[]
    for keys, sub in df.groupby(list(group_cols)):
        if len(sub) >= 4:
            stat,p = shapiro(sub['EditingValue'])
            rows.append({**{c:k for c,k in zip(group_cols,keys)}, 'n':len(sub), 'shapiro_p':p})
    return pd.DataFrame(rows)

def mannwhitney_ex_vs_sed(df):
    """Runs at two levels: by Phase including all ZT, and by Phase×ZT.
       Returns a DataFrame with p and p_FDR."""
    out=[]

    # 1) By Phase (all ZTs combined)
    for phase, sub in df.groupby('Phase'):
        a = sub[sub['treatment']=='Exercise']['EditingValue']
        b = sub[sub['treatment']=='Sedentary']['EditingValue']
        if len(a)>=3 and len(b)>=3:
            stat,p = mannwhitneyu(a,b,alternative='two-sided')
            out.append({'level':'Phase_only','Phase':phase,'zt':None,'n_ex':len(a),'n_sed':len(b),'U':stat,'p':p})

    # 2) By Phase × ZT
    for (phase, zt), sub in df.groupby(['Phase','zt']):
        a = sub[sub['treatment']=='Exercise']['EditingValue']
        b = sub[sub['treatment']=='Sedentary']['EditingValue']
        if len(a)>=3 and len(b)>=3:
            stat,p = mannwhitneyu(a,b,alternative='two-sided')
            out.append({'level':'Phase_ZT','Phase':phase,'zt':zt,'n_ex':len(a),'n_sed':len(b),'U':stat,'p':p})

    res = pd.DataFrame(out).sort_values(['level','Phase','zt'], ignore_index=True)
    if not res.empty:
        # FDR correction for all tests in this file
        res['p_fdr'] = multipletests(res['p'], method='fdr_bh')[1]
        # Map to significance stars
        def stars(p):
            return 'ns' if p>=0.05 else ('*' if p<0.05 and p>=0.01 else ('**' if p<0.01 and p>=0.001 else '***'))
        res['sig'] = res['p_fdr'].apply(stars)
    return res

def annotate_sig_on_axis(ax, xpos, top, text):
    """Add significance stars above the box at xpos. 'top' is the height."""
    ax.text(xpos, top, text, ha='center', va='bottom', fontsize=11)

def add_sig_to_phase_plot(ax, df, res_phase_only, x_order):
    """Plot 1: Phase x Treatment. Adds significance star for each Phase comparing Exercise vs Sedentary."""
    y = df['EditingValue']
    y_min, y_max = np.min(y), np.max(y)
    pad = (y_max - y_min) * 0.08 if y_max>y_min else 0.02

    for i, phase in enumerate(x_order):
        # Set height above the highest box for this Phase
        ymax_phase = df[df['Phase']==phase]['EditingValue'].max()
        top = ymax_phase + pad
        row = res_phase_only[(res_phase_only['level']=='Phase_only') & (res_phase_only['Phase']==phase)]
        if not row.empty:
            annotate_sig_on_axis(ax, i, top, row['sig'].iloc[0])

def add_sig_to_zt_treat_plot(ax, df, res_phase_zt, x_order):
    """Plot 3: ZT x Treatment. Adds significance star for each ZT (aggregated across Phases).
       If you prefer by each Phase separately – use Plot 4."""
    y = df['EditingValue']
    y_min, y_max = np.min(y), np.max(y)
    pad = (y_max - y_min) * 0.08 if y_max>y_min else 0.02

    # Combine both Phases for each ZT: recalculate test
    tmp=[]
    for zt, sub in df.groupby('zt'):
        a = sub[sub['treatment']=='Exercise']['EditingValue']
        b = sub[sub['treatment']=='Sedentary']['EditingValue']
        if len(a)>=3 and len(b)>=3:
            stat,p = mannwhitneyu(a,b,alternative='two-sided')
            tmp.append({'zt':zt,'p':p})
    if tmp:
        tmp = pd.DataFrame(tmp).sort_values('zt')
        tmp['p_fdr'] = multipletests(tmp['p'], method='fdr_bh')[1]
        def stars(p): return 'ns' if p>=0.05 else ('*' if p<0.05 and p>=0.01 else ('**' if p<0.01 and p>=0.001 else '***'))
        tmp['sig'] = tmp['p_fdr'].apply(stars)

        for i, zt in enumerate(x_order):
            ymax_zt = df[df['zt']==zt]['EditingValue'].max()
            top = ymax_zt + pad
            row = tmp[tmp['zt']==zt]
            if not row.empty:
                annotate_sig_on_axis(ax, i, top, row['sig'].iloc[0])

def add_sig_to_facetgrid(g, df, res_phase_zt):
    """Plot 4: Facet by Phase. Adds significance stars for each ZT within the Phase."""
    for ax, (phase, subdf) in zip(g.axes.flat, df.groupby('Phase')):
        y = subdf['EditingValue']
        pad = (y.max()-y.min())*0.08 if y.max()>y.min() else 0.02
        x_order = sorted(subdf['zt'].unique())
        for i, zt in enumerate(x_order):
            ymax_zt = subdf[subdf['zt']==zt]['EditingValue'].max()
            top = ymax_zt + pad
            row = res_phase_zt[(res_phase_zt['level']=='Phase_ZT') & (res_phase_zt['Phase']==phase) & (res_phase_zt['zt']==zt)]
            if not row.empty:
                annotate_sig_on_axis(ax, i, top, row['sig'].iloc[0])

def run_all_for_file(csv_filename):
    file_label = csv_filename.split("_")[0]
    df = pd.read_csv(csv_filename)
    df = rename_edit_col(df)
    base_filename = os.path.splitext(os.path.basename(csv_filename))[0]

    # ---------- Statistics ----------
    shapiro_df = shapiro_by_group(df)
    mw_df = mannwhitney_ex_vs_sed(df)

    # Save result files
    shapiro_path = f"{base_filename}_shapiro_by_group.csv"
    mw_path = f"{base_filename}_mannwhitney_ex_vs_sed_FDR.csv"
    shapiro_df.to_csv(shapiro_path, index=False)
    mw_df.to_csv(mw_path, index=False)

    print(f"Saved:\n  {shapiro_path}\n  {mw_path}")

    # ---------- Plots ----------
    sns.set(style="whitegrid")

    # Plot 1: Phase on X-axis, hue by Treatment
    plt.figure(figsize=(8, 6))
    ax = sns.boxplot(data=df, x="Phase", y="EditingValue", hue="treatment")
    sns.stripplot(data=df, x="Phase", y="EditingValue", hue="treatment",
                  dodge=True, color='black', alpha=0.3, jitter=True)
    # Ensure single legend
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        plt.legend(handles[:2], ['Exercise','Sedentary'], title="Treatment", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.suptitle("A2G Editing Index by Phase and Treatment", fontsize=14)
    plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')
    plt.ylabel("A2G Editing Index"); plt.xlabel("Phase")
    plt.tight_layout()

    # Add significance stars by Phase only
    x_order = list(df['Phase'].drop_duplicates())
    add_sig_to_phase_plot(ax, df, mw_df, x_order)

    out1 = f"{base_filename}_A2G_by_phase_and_treatment_SIG.png"
    plt.savefig(out1, dpi=300); plt.close()

    # Plot 2: ZT on X-axis, hue by Phase (no Treatment comparison, leave without stars)
    plt.figure(figsize=(12, 6))
    ax = sns.boxplot(data=df, x="zt", y="EditingValue", hue="Phase")
    sns.stripplot(data=df, x="zt", y="EditingValue", hue="Phase",
                  dodge=True, color='black', alpha=0.2, jitter=True)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        plt.legend(handles[:2], ['Rest phase (ZT3)','Active phase (ZT15)'], title="Phase",
                   bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.suptitle("A2G Editing Index by ZT and Phase (All Treatments)", fontsize=14)
    plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')
    plt.ylabel("A2G Editing Index"); plt.xlabel("Zeitgeber Time (ZT)")
    plt.tight_layout()
    out2 = f"{base_filename}_A2G_all_by_zt_phase.png"
    plt.savefig(out2, dpi=300); plt.close()

    # Plot 3: ZT on X-axis, hue by Treatment – add stars for each ZT (aggregated Phases)
    plt.figure(figsize=(12, 6))
    ax = sns.boxplot(data=df, x="zt", y="EditingValue", hue="treatment")
    sns.stripplot(data=df, x="zt", y="EditingValue", hue="treatment",
                  dodge=True, color='black', alpha=0.3, jitter=True)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        plt.legend(handles[:2], ['Exercise','Sedentary'], title="Treatment", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.suptitle("A2G Editing Index by ZT and Treatment", fontsize=14)
    plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')
    plt.ylabel("A2G Editing Index"); plt.xlabel("Zeitgeber Time (ZT)")
    plt.tight_layout()

    x_order = sorted(df['zt'].unique())
    add_sig_to_zt_treat_plot(ax, df, mw_df, x_order)

    out3 = f"{base_filename}_A2G_by_zt_and_treatment_SIG.png"
    plt.savefig(out3, dpi=300); plt.close()

    # Plot 4: Facet by Phase, hue by Treatment – stars for each ZT inside each Phase
    g = sns.catplot(
        data=df, x="zt", y="EditingValue",
        hue="treatment", col="Phase",
        kind="box", height=5, aspect=1.2
    )
    # Add points
    for ax, (phase, subdf) in zip(g.axes.flat, df.groupby('Phase')):
        sns.stripplot(data=subdf, x='zt', y='EditingValue', hue='treatment', dodge=True,
                      color='black', alpha=0.3, jitter=True, ax=ax, legend=False)
    g.fig.subplots_adjust(top=0.85)
    g.figure.suptitle("A2G Editing Index by ZT and Treatment – Split by Phase", fontsize=14)
    g.figure.text(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')

    # Ensure single legend only
    g.add_legend(title='Treatment')

    add_sig_to_facetgrid(g, df, mw_df)

    out4 = f"{base_filename}_A2G_by_zt_and_treatment_split_by_phase_SIG.png"
    g.savefig(out4, dpi=300); plt.close(g.figure)

    return {
        'shapiro_csv': shapiro_path,
        'mannwhitney_csv': mw_path,
        'figs':[out1,out2,out3,out4]
    }

# -----------------------
# Run on both datasets
# -----------------------
results_alu = run_all_for_file("Alu_EditingIndex_with_groups_A2G_only.csv")
results_utr3 = run_all_for_file("3UTR_EditingIndex_with_groups_A2G_only.csv")

print("ALU:", results_alu)
print("3UTR:", results_utr3)
