"""
扇环形工具模块 (annular_sector)

本模块提供对“偏移扇环形”的几何计算，所有公开函数均用于计算将 n 等分扇环形
的两条径向边界向内水平偏移距离 b 后所得曲边四边形的面积和周长。
内部辅助函数遵循单一职责原则，供公开函数调用。

核心概念：
- 扇环形：由两个同心圆（内半径 r，外半径 R，0<r<R）和圆心角 a（度）围成。
- 偏移扇环形：将上述扇环形的两条径向边界同时向内水平移动距离 b（0<b<r）
  后，与外圆弧、内圆弧及两条偏移边（直线段）围成的曲边四边形。
- **旋转处理**：为简化一侧的计算，将整个扇环形旋转，使其一条原始径向边与 y 轴重合。
  这样，向内水平偏移 b 后，该侧边界变为直线 x = b，它与圆的交点位于第一象限，
  便于用直角坐标解析。另一侧由于对称性（偏移量相同）会产生完全相同的几何量，
  因此只需计算一侧，然后加倍即可。
- **截止角 (cutoff angle)**：θ(ρ) = arcsin(b/ρ)（弧度），表示在上述旋转后的坐标系中，
  直线 x = b 与半径为 ρ 的圆的交点与原点连线相对于 y 轴的夹角。该角即为偏移过程中
  该侧径向边上的点沿圆弧移动所扫过的圆心角。由于对称性，另一侧同样产生 θ(ρ)。
- 原始圆心角：a = 360/n（度）。
- 剩余圆心角：a_ρ' = a - 2·(180/π)·θ(ρ)（度），即偏移后半径为 ρ 的圆弧
  对应的圆心角（由原始角减去两个对称的截止角得到）。
- 原始扇环面积：S = (a/360) * π * (R² - r²)。
- 切去扇环面积：旋转后，一侧切去的区域是由 x=0 到 x=b 的两条竖直线在同心圆间
  截出的扇环形，其面积记为 S_cut = _area_of_offset(r, R, b)。由于对称性，
  另一侧也切去同样大小的面积。
- 偏移扇环面积：S_offset = S - 2 * S_cut。
- 偏移边边长：L = √(R² - b²) - √(r² - b²)，由旋转后坐标系中直线 x = b 被两圆截得的线段长度给出。
- 偏移扇环周长：P_offset = 2 * L + 内圆弧长 + 外圆弧长（两倍的直线段加上对称的两段圆弧）。

所有函数均使用纯 Python 和 math 库实现，避免第三方依赖。
"""

import math

# ---------- 私有辅助函数 ----------


def _annular_sector_area(inner: float, outer: float, angle_deg: float) -> float:
    """
    计算原始扇环形面积。

    参数:
        inner:     内圆半径，必须小于 outer
        outer:     外圆半径，必须大于 inner
        angle_deg: 圆心角，以度为单位，必须在 (0, 360] 范围内

    返回:
        扇环形面积 = (angle_deg/360) * π * (outer² - inner²)
    """
    if not (0 < inner < outer):
        raise ValueError("半径必须满足 0 < inner < outer")
    if not (0 < angle_deg <= 360):
        raise ValueError("圆心角 angle_deg 应在 (0, 360] 度之间")
    return (angle_deg / 360.0) * math.pi * (outer * outer - inner * inner)


def _curved_trapezoid_area(radius: float, offset: float) -> float:
    """
    计算由圆 x² + y² = radius²、直线 y = 0、x = 0、x = offset
    在第一象限内围成的曲边梯形面积。

    参数:
        radius: 圆的半径，必须满足 0 < offset < radius
        offset: x 坐标右边界（水平偏移量）

    返回:
        曲边梯形面积
    """
    if not (0 < offset < radius):
        raise ValueError("参数必须满足 0 < offset < radius")
    term1 = (offset / 2.0) * ((radius * radius - offset * offset) ** 0.5)
    term2 = (radius * radius / 2.0) * math.asin(offset / radius)
    return term1 + term2


def _area_of_offset(inner: float, outer: float, offset: float) -> float:
    """
    计算由两个同心圆 (x²+y²=inner² 和 x²+y²=outer²) 与两条竖直线
    (x=0, x=offset) 在第一象限内围成的扇环形面积。
    这是在旋转后坐标系中，单侧切去的一块面积。

    参数:
        inner:  内圆半径，必须满足 0 < offset < inner < outer
        outer:  外圆半径，必须大于 inner
        offset: 水平偏移量

    返回:
        单侧切去的扇环形面积
    """
    if not (0 < offset < inner < outer):
        raise ValueError("参数必须满足 0 < offset < inner < outer")
    return _curved_trapezoid_area(outer, offset) - _curved_trapezoid_area(inner, offset)


def _offset_edge_length(inner: float, outer: float, offset: float) -> float:
    """
    计算直线 x = offset 在第一象限内被两个同心圆 x²+y²=inner² 和 x²+y²=outer²
    所截的线段长度。这是在旋转后坐标系中，单条偏移边的边长。

    参数:
        inner:  内圆半径，必须满足 0 < offset < inner < outer
        outer:  外圆半径，必须大于 inner
        offset: 直线 x = offset 的 x 坐标

    返回:
        单条偏移边边长 L = √(outer² - offset²) - √(inner² - offset²)
    """
    if not (0 < offset < inner < outer):
        raise ValueError("参数必须满足 0 < offset < inner < outer")
    return (outer * outer - offset * offset) ** 0.5 - (
        inner * inner - offset * offset
    ) ** 0.5


def _cutoff_angle(radius: float, offset: float) -> float:
    """
    计算截止角（弧度）。此计算基于旋转后的坐标系，其中一条原始径向边与 y 轴重合。
    偏移后得到直线 x = offset，它与圆 x²+y²=radius² 的交点位于第一象限。
    该交点与原点的连线与 y 轴的夹角即为截止角，它等于偏移过程中该侧径向边上的点
    沿圆弧扫过的圆心角。几何关系：sinθ = offset/radius。

    参数:
        radius: 圆的半径，必须满足 0 < offset < radius
        offset: 水平偏移量

    返回:
        截止角 θ，单位为弧度
    """
    if not (0 < offset < radius):
        raise ValueError("参数必须满足 0 < offset < radius")
    return math.asin(offset / radius)


def _remaining_angle_after_offset(
    radius: float, angle_deg: float, offset: float
) -> float:
    """
    计算将圆心角为 angle_deg 的扇形（半径为 radius）的两条径向边界向内水平偏移 offset 后，
    在半径为 radius 的圆弧上所截得的圆弧对应的圆心角（度）。

    由于两条径向边对称偏移，每条边产生的截止角均为 θ = _cutoff_angle(radius, offset)（弧度）。
    因此，剩余圆心角 = angle_deg - 2·(180/π)·θ（度）。

    参数:
        radius:    圆的半径，必须满足 0 < offset < radius
        angle_deg: 原始扇形的圆心角，以度为单位
        offset:    水平偏移量

    返回:
        剩余圆心角，以度为单位

    异常:
        ValueError: 如果剩余圆心角 ≤ 0，表示偏移量过大或原始角过小
    """
    if not (0 < offset < radius):
        raise ValueError("参数必须满足 0 < offset < radius")
    theta_rad = _cutoff_angle(radius, offset)
    theta_deg = math.degrees(theta_rad)
    remaining_angle_deg = angle_deg - 2.0 * theta_deg
    if remaining_angle_deg <= 0:
        raise ValueError("偏移量过大或原始角度过小，导致剩余圆弧非正")
    return remaining_angle_deg


def _arc_length_after_offset(radius: float, angle_deg: float, offset: float) -> float:
    """
    计算将圆心角为 angle_deg 的扇形（半径为 radius）的两条径向边界向内水平偏移 offset 后，
    在半径为 radius 的圆弧上所截得的弧长。

    几何意义：先由 _remaining_angle_after_offset 计算剩余圆心角，再乘以半径得到弧长。

    参数:
        radius:    圆的半径，必须满足 0 < offset < radius
        angle_deg: 原始扇形的圆心角，以度为单位
        offset:    水平偏移量

    返回:
        偏移后的圆弧长
    """
    remaining_angle_deg = _remaining_angle_after_offset(radius, angle_deg, offset)
    return math.radians(remaining_angle_deg) * radius


# ---------- 公开函数 ----------


def offset_area(r: float, R: float, n: int, b: float) -> float:
    """
    计算将 n 等分扇环形（圆心角 = 360/n 度）的两条径向边界同时向内
    水平偏移距离 b 后所得曲边四边形的面积。

    面积由原始扇环面积减去两个对称切去的扇环面积（每个由旋转后坐标系中
    x=0 到 x=b 围成）得到。

    参数:
        r: 内圆半径，必须满足 0 < b < r < R（但内部会自动排序）
        R: 外圆半径，必须大于 r（内部会自动交换以保证 inner < outer）
        n: 等分份数，必须为正整数（例如 n=8 表示八等分）
        b: 水平偏移量，必须满足 0 < b < min(r,R)

    返回:
        偏移扇环形的面积
    """
    # 自动排序半径
    inner = min(r, R)
    outer = max(r, R)
    if not (0 < b < inner):
        raise ValueError("参数必须满足 0 < b < 内圆半径")
    if n <= 0:
        raise ValueError("等分份数 n 必须为正数")

    angle_deg = 360.0 / n
    original_area = _annular_sector_area(inner, outer, angle_deg)
    cut_area = _area_of_offset(inner, outer, b)
    return original_area - 2.0 * cut_area


def offset_perimeter(r: float, R: float, n: int, b: float) -> float:
    """
    计算将 n 等分扇环形（圆心角 = 360/n 度）的两条径向边界同时向内
    水平偏移距离 b 后所得曲边四边形的周长。

    周长由两倍的偏移边边长（对称）加上内外两段圆弧长组成。

    参数:
        r: 内圆半径，必须满足 0 < b < r < R（内部自动排序）
        R: 外圆半径，必须大于 r
        n: 等分份数，必须为正整数
        b: 水平偏移量，必须满足 0 < b < min(r,R)

    返回:
        偏移扇环形的周长
    """
    inner = min(r, R)
    outer = max(r, R)
    if not (0 < b < inner):
        raise ValueError("参数必须满足 0 < b < 内圆半径")
    if n <= 0:
        raise ValueError("等分份数 n 必须为正数")

    angle_deg = 360.0 / n

    # 检查内圆上是否有正的剩余圆弧（最严格的条件）
    # 利用 _remaining_angle_after_offset 检查，但此处只需验证内圆即可
    # 直接调用 _remaining_angle_after_offset 会抛出异常，但我们需要先验证有效性
    # 为了一致性，我们仍使用 _cutoff_angle 检查，避免重复计算
    inner_cutoff_rad = _cutoff_angle(inner, b)
    if angle_deg <= 2.0 * math.degrees(inner_cutoff_rad):
        raise ValueError("偏移量过大或等分数过小，导致内圆上无剩余圆弧")

    edge_len = _offset_edge_length(inner, outer, b)
    inner_arc = _arc_length_after_offset(inner, angle_deg, b)
    outer_arc = _arc_length_after_offset(outer, angle_deg, b)

    return 2.0 * edge_len + inner_arc + outer_arc


# ---------- 示例测试 ----------
if __name__ == "__main__":
    r_val, R_val, n_val, b_val = 3270, 3475, 16, 100
    area = offset_area(r_val, R_val, n_val, b_val)
    perimeter = offset_perimeter(r_val, R_val, n_val, b_val)
    print(f"偏移扇环形面积   = {area}")
    print(f"偏移扇环形周长   = {perimeter}")
