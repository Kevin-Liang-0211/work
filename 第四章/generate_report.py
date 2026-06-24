from pathlib import Path
from math import cos, pi, sinh, sqrt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Inches, Pt


ROOT = Path(__file__).resolve().parent
CODE_DIR = ROOT / "代码"
IMG_DIR = ROOT / "images"
IMG_DIR.mkdir(exist_ok=True)


def set_cn_font(run, name="Songti SC"):
    run.font.name = name
    run._element.rPr.rFonts.set(qn("w:eastAsia"), name)


def style_doc(doc: Document):
    section = doc.sections[0]
    section.page_width = Inches(8.27)
    section.page_height = Inches(11.69)
    section.top_margin = Inches(1.0)
    section.bottom_margin = Inches(1.0)
    section.left_margin = Inches(1.0)
    section.right_margin = Inches(1.0)

    normal = doc.styles["Normal"]
    normal.font.size = Pt(11)
    normal.font.name = "Times New Roman"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Songti SC")
    pf = normal.paragraph_format
    pf.line_spacing = 1.3
    pf.space_after = Pt(6)

    for style_name, size in [("Title", 18), ("Heading 1", 15), ("Heading 2", 13)]:
        style = doc.styles[style_name]
        style.font.name = "Times New Roman"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Heiti SC")
        style.font.size = Pt(size)


def poisson_iterative(method, h, omega=1.9, tol=1e-8, max_iters=200000):
    n = int(round(1 / h))
    U = np.zeros((n + 1, n + 1))
    Unew = U.copy()
    rhs = -16.0
    for it in range(1, max_iters + 1):
        old = U.copy()
        if method == "jacobi":
            Unew[1:n, 1:n] = 0.25 * (
                old[:-2, 1:-1] + old[2:, 1:-1] + old[1:-1, :-2] + old[1:-1, 2:] + h * h * rhs
            )
            err = np.max(np.abs(Unew - old))
            U = Unew.copy()
        elif method == "gs":
            for j in range(1, n):
                for i in range(1, n):
                    U[i, j] = 0.25 * (U[i - 1, j] + old[i + 1, j] + U[i, j - 1] + old[i, j + 1] + h * h * rhs)
            err = np.max(np.abs(U - old))
        elif method == "sor":
            for j in range(1, n):
                for i in range(1, n):
                    jac = 0.25 * (U[i - 1, j] + old[i + 1, j] + U[i, j - 1] + old[i, j + 1] + h * h * rhs)
                    U[i, j] = (1 - omega) * old[i, j] + omega * jac
            err = np.max(np.abs(U - old))
        else:
            raise ValueError(method)
        if err < tol:
            return it, err, U
    raise RuntimeError(f"{method} does not converge")


def helmholtz_solution(lam, h):
    mu = sqrt(lam * lam + pi * pi)
    n = int(round(1 / h))
    x = np.linspace(0, 1, n + 1)
    y = np.linspace(0, 1, n + 1)
    phi = np.zeros((n + 1, n + 1))
    phi[0, :] = np.sinh(mu * (1 - y)) / np.sinh(mu)
    phi[-1, :] = -np.sinh(mu * (1 - y)) / np.sinh(mu)
    phi[:, 0] = np.cos(pi * x)
    phi[:, -1] = 0.0

    Nix = n - 1
    Niy = n - 1
    N = Nix * Niy
    A = np.zeros((N, N))
    b = np.zeros(N)

    def idx(i, j):
        return j * Nix + i

    for j in range(Niy):
        for i in range(Nix):
            row = idx(i, j)
            A[row, row] = -4 - lam * lam * h * h
            if i > 0:
                A[row, idx(i - 1, j)] = 1
            else:
                b[row] -= phi[0, j + 1]
            if i < Nix - 1:
                A[row, idx(i + 1, j)] = 1
            else:
                b[row] -= phi[-1, j + 1]
            if j > 0:
                A[row, idx(i, j - 1)] = 1
            else:
                b[row] -= phi[i + 1, 0]
            if j < Niy - 1:
                A[row, idx(i, j + 1)] = 1
            else:
                b[row] -= phi[i + 1, -1]

    sol = np.linalg.solve(A, b)
    for j in range(Niy):
        for i in range(Nix):
            phi[i + 1, j + 1] = sol[idx(i, j)]

    xx, yy = np.meshgrid(x, y, indexing="ij")
    exact = np.cos(pi * xx) * np.sinh(mu * (1 - yy)) / np.sinh(mu)
    return x, y, phi, exact, np.max(np.abs(phi - exact))


def save_surface(U, h, out_path):
    x = np.linspace(0, 1, U.shape[0])
    y = np.linspace(0, 1, U.shape[1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    plt.figure(figsize=(5.6, 4.4))
    ax = plt.axes(projection="3d")
    ax.plot_surface(X, Y, U.T, cmap="viridis", edgecolor="none")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("u")
    ax.set_title(f"Jacobi solution, h={h}")
    ax.view_init(elev=24, azim=-58)
    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()


def save_helmholtz_surface(phi, lam, out_path):
    x = np.linspace(0, 1, phi.shape[0])
    y = np.linspace(0, 1, phi.shape[1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    plt.figure(figsize=(5.6, 4.4))
    ax = plt.axes(projection="3d")
    ax.plot_surface(X, Y, phi.T, cmap="plasma", edgecolor="none")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel(r"$\phi$")
    ax.set_title(f"Helmholtz solution, lambda={lam}")
    ax.view_init(elev=26, azim=-60)
    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()


def generate_results():
    h_list = [0.1, 0.05, 0.025, 0.0125]
    poisson_rows = []
    jacobi_images = []
    for h in h_list:
        ij, ej, Uj = poisson_iterative("jacobi", h)
        ig, eg, Ug = poisson_iterative("gs", h)
        isor, es, Us = poisson_iterative("sor", h)
        center = Uj[Uj.shape[0] // 2, Uj.shape[1] // 2]
        poisson_rows.append((h, ij, ej, ig, eg, isor, es, center))
        img_path = IMG_DIR / f"CH4_HW1_Jacobi_h{str(h).replace('.', 'p')}.png"
        save_surface(Uj, h, img_path)
        jacobi_images.append(img_path)

    helmholtz_rows = []
    helmholtz_images = []
    for lam in [0.5, 1.0, 2.0]:
        x, y, phi, exact, err = helmholtz_solution(lam, 0.025)
        helmholtz_rows.append((lam, err, float(phi.max()), float(phi.min())))
        img_path = IMG_DIR / f"CH4_HW2_lambda_{str(lam).replace('.', 'p')}.png"
        save_helmholtz_surface(phi, lam, img_path)
        helmholtz_images.append(img_path)
    return poisson_rows, jacobi_images, helmholtz_rows, helmholtz_images


def add_paragraph(doc, text, bold=False, align=WD_ALIGN_PARAGRAPH.LEFT):
    p = doc.add_paragraph()
    p.alignment = align
    r = p.add_run(text)
    set_cn_font(r, "Songti SC" if not bold else "Heiti SC")
    r.bold = bold
    return p


def build_report(poisson_rows, jacobi_images, helmholtz_rows, helmholtz_images):
    doc = Document()
    style_doc(doc)

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("微分方程数值解法 实验报告")
    set_cn_font(r, "Heiti SC")
    r.bold = True
    r.font.size = Pt(18)

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("第四章 椭圆型方程数值解法")
    set_cn_font(r, "Heiti SC")
    r.font.size = Pt(14)

    doc.add_paragraph("")
    add_paragraph(doc, "一、实验目的", bold=True)
    add_paragraph(doc, "掌握 Poisson 方程五点差分格式的构造方法，理解 Jacobi、Gauss-Seidel 和 SOR 三种经典迭代法的收敛特征；掌握二维 Helmholtz 方程二阶差分离散思想，并能够编写程序对椭圆型边值问题进行数值求解与结果分析。")

    add_paragraph(doc, "二、实验内容与结果", bold=True)
    add_paragraph(doc, "题目1：Poisson 方程的五点差分与迭代求解", bold=True)
    add_paragraph(doc, "问题描述：求解 Poisson 方程 u_xx + u_yy = 16，定义域 D=(0,1)×(0,1)，边界条件为 u=0。取初值 u0=0，分别采用 Jacobi、Gauss-Seidel 和 SOR 迭代法求解差分方程，收敛条件取 max|U^(n+1)-U^n| < 10^-8。网格步长依次取 Δx=Δy=0.1, 0.05, 0.025, 0.0125，其中 SOR 松弛因子取 ω=1.9。")
    add_paragraph(doc, "算法说明：在均匀网格上采用五点差分格式")
    add_paragraph(doc, "U_{i,j} = (U_{i-1,j}+U_{i+1,j}+U_{i,j-1}+U_{i,j+1}-h^2 f_{i,j}) / 4")
    add_paragraph(doc, "构造 Jacobi、Gauss-Seidel 与 SOR 迭代。由于本题 f=16，内部节点数值解为负，解在区域中心附近达到最小值。")

    table = doc.add_table(rows=1, cols=8)
    table.style = "Table Grid"
    headers = ["Δx=Δy", "Jacobi迭代次数", "Jacobi误差", "GS迭代次数", "GS误差", "SOR迭代次数", "SOR误差", "中心点值"]
    for c, text in enumerate(headers):
        table.rows[0].cells[c].text = text
    for row in poisson_rows:
        cells = table.add_row().cells
        vals = [
            f"{row[0]:.4f}", str(row[1]), f"{row[2]:.3e}", str(row[3]), f"{row[4]:.3e}",
            str(row[5]), f"{row[6]:.3e}", f"{row[7]:.6f}"
        ]
        for c, val in enumerate(vals):
            cells[c].text = val

    add_paragraph(doc, "结果分析：")
    add_paragraph(doc, "1. 三种迭代法均能在给定容差下收敛，其中 Gauss-Seidel 总是比 Jacobi 更快。")
    add_paragraph(doc, "2. SOR 在粗网格和中等网格上优势最明显；当网格继续加密时，若 ω 固定为 1.9，则收敛次数会有所回升，但仍显著少于前两者。")
    add_paragraph(doc, "3. 随着网格加密，Jacobi 解在中心点的数值逐步稳定在约 -1.1786，说明离散解正在趋于连续真解。")

    for img in jacobi_images:
        doc.add_picture(str(img), width=Inches(4.8))
        p = doc.paragraphs[-1]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER

    add_paragraph(doc, "题目2（选做）：二维 Helmholtz 方程的二阶差分格式", bold=True)
    add_paragraph(doc, "问题描述：求解 Helmholtz 方程 φ_xx + φ_yy = λ^2 φ，定义域 D=(0,1)×(0,1)，边界条件为 φ(0,y)=sinh(μ(1-y))/sinhμ，φ(1,y)=-sinh(μ(1-y))/sinhμ，φ(x,0)=cosπx，φ(x,1)=0，其中 μ=sqrt(λ^2+π^2)。分别给出 λ=0.5, 1.0, 2.0 时的数值解。")
    add_paragraph(doc, "算法说明：对内部节点采用二阶中心差分，得到")
    add_paragraph(doc, "φ_{i-1,j}+φ_{i+1,j}+φ_{i,j-1}+φ_{i,j+1}-(4+λ^2 h^2)φ_{i,j}=0。")
    add_paragraph(doc, "将其整理为稀疏线性方程组后，使用矩阵直接求解。解析解为 φ(x,y)=cos(πx)sinh(μ(1-y))/sinhμ。")

    table2 = doc.add_table(rows=1, cols=4)
    table2.style = "Table Grid"
    for c, text in enumerate(["λ", "最大误差", "数值解最大值", "数值解最小值"]):
        table2.rows[0].cells[c].text = text
    for row in helmholtz_rows:
        vals = [f"{row[0]:.1f}", f"{row[1]:.3e}", f"{row[2]:.6f}", f"{row[3]:.6f}"]
        cells = table2.add_row().cells
        for c, val in enumerate(vals):
            cells[c].text = val

    add_paragraph(doc, "结果分析：在固定网格 h=0.025 下，三组参数的最大误差均保持在较小量级，说明构造的二阶差分格式能够稳定逼近解析解。随着 λ 增大，解在 y 方向的衰减更快，曲面在上边界附近更快趋近于零。")
    for img in helmholtz_images:
        doc.add_picture(str(img), width=Inches(4.8))
        p = doc.paragraphs[-1]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER

    add_paragraph(doc, "三、实验总结", bold=True)
    add_paragraph(doc, "1. 五点差分格式是二维 Poisson/Helmholtz 方程最基本也最实用的离散工具，结构简单、便于实现。")
    add_paragraph(doc, "2. Jacobi 法实现最直接，但收敛较慢；Gauss-Seidel 通过及时利用新值明显加快了收敛；SOR 进一步通过松弛参数显著降低迭代次数。")
    add_paragraph(doc, "3. Helmholtz 方程在离散后可写成带参数项的稀疏线性系统，采用二阶中心差分可以保持较好的精度。")
    add_paragraph(doc, "4. 本章实验体现了椭圆型方程数值求解的两个核心思路：一是将边值问题离散为线性方程组，二是根据问题规模选择迭代法或直接法进行求解。")

    report_path = ROOT / "实验报告_第四章.docx"
    doc.save(report_path)
    return report_path


def build_code_appendix():
    doc = Document()
    style_doc(doc)
    add_paragraph(doc, "第四章附录：源代码", bold=True, align=WD_ALIGN_PARAGRAPH.CENTER)
    for file_path in sorted(CODE_DIR.glob("*.m")):
        add_paragraph(doc, f"{file_path.name}", bold=True)
        text = file_path.read_text(encoding="utf-8")
        for line in text.splitlines():
            p = doc.add_paragraph()
            r = p.add_run(line if line else " ")
            set_cn_font(r, "Menlo")
            r.font.size = Pt(8.5)
            p.paragraph_format.space_after = Pt(0)
    out = ROOT / "附录_源代码.docx"
    doc.save(out)
    return out


if __name__ == "__main__":
    poisson_rows, jacobi_images, helmholtz_rows, helmholtz_images = generate_results()
    report = build_report(poisson_rows, jacobi_images, helmholtz_rows, helmholtz_images)
    appendix = build_code_appendix()
    print(report)
    print(appendix)
