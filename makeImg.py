
import matplotlib.pyplot as plt
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def plot_plddt(cif_path, chain_id=None, highlight_ranges=None):
    """
    AlphaFold 3 CIF 파일에서 pLDDT를 추출하고 그래프로 출력한다.

    highlight_ranges:
        [(17, 19), (32, 34), (41, 42)]
    """

    cif = MMCIF2Dict(cif_path)

    atom_names = cif["_atom_site.label_atom_id"]
    chain_ids = cif["_atom_site.label_asym_id"]
    residue_ids = cif["_atom_site.label_seq_id"]
    residue_names = cif["_atom_site.label_comp_id"]
    b_factors = cif["_atom_site.B_iso_or_equiv"]

    data = []

    for atom, chain, residue, res_name, b_factor in zip(
        atom_names,
        chain_ids,
        residue_ids,
        residue_names,
        b_factors
    ):
        if atom != "CA":
            continue

        if chain_id is not None and chain != chain_id:
            continue

        data.append({
            "Chain": chain,
            "Residue": int(residue),
            "Residue_Name": res_name,
            "pLDDT": float(b_factor)
        })

    if not data:
        raise ValueError("C-alpha 원자를 찾을 수 없습니다.")

    df = pd.DataFrame(data)

    for chain, group in df.groupby("Chain"):

        group = group.sort_values("Residue")

        
        
        # 그래프 생성
        fig, ax = plt.subplots(figsize=(16, 8))

        # 색상 설정
        line_color = "#1f77b4"
        highlight_color = "#1f77b4"

        # 흰색 배경
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        # pLDDT 그래프
        ax.plot(
            group["Residue"],
            group["pLDDT"],
            color=line_color,
            linewidth=1.8,
            marker="o",
            markersize=5
        )

        # pLDDT 기준선
        for value in [50, 70, 90]:
            ax.axhline(
                y=value,
                color=line_color,
                linestyle="--",
                linewidth=1.5,
                alpha=0.9
            )

        # 구간 강조
        if highlight_ranges:
            for start, end in highlight_ranges:

                ax.axvspan(
                    start - 0.5,
                    end + 0.5,
                    color=highlight_color,
                    alpha=0.12,
                    zorder=0
                )

                # 구간 이름 표시
                ax.text(
                    (start + end) / 2,
                    93,
                    f"{start}-{end}",
                    ha="center",
                    va="center",
                    fontsize=12,
                    color="black",
                    clip_on=True
                )

        # 축 범위
        ax.set_ylim(0, 100)

        ax.set_xlim(
            group["Residue"].min() - 0.5,
            group["Residue"].max() + 0.5
        )

        # 축 제목
        ax.set_xlabel(
            "α-synuclein residue position",
            fontsize=13
        )

        ax.set_ylabel(
            "Cα pLDDT",
            fontsize=13
        )

        # 그래프 제목
        ax.set_title(
            f"AlphaFold Server prediction: α-synuclein residue-level pLDDT (Model 0)",
            fontsize=16,
            pad=10
        )

        # 격자
        ax.grid(
            True,
            color="#cccccc",
            alpha=0.4,
            linewidth=0.8
        )

        # 눈금
        ax.tick_params(
            axis="both",
            labelsize=11
        )

        # 테두리
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(1.0)

        plt.tight_layout()
        plt.show()



    return df


# 실행
df = plot_plddt(
    "model0.cif",
    highlight_ranges=[
        (1, 60),
        (61, 95),
        (96, 140)
    ]
)

df.to_csv("plddt_results.csv", index=False)
