"""
算术平均和对数平均计算模块 (util.mean)

提供算术平均、对数平均以及它们的几何平均的计算函数。
所有函数均使用纯 Python 和 math 库实现，避免第三方依赖。
"""

import math


def arithmetic_mean(a: float, b: float) -> float:
    """
    计算两个数的算术平均值。

    算术平均定义为 (a + b) / 2。

    参数:
        a: 第一个数
        b: 第二个数

    返回:
        算术平均值
    """
    return (a + b) / 2.0


def logarithmic_mean(a: float, b: float) -> float:
    """
    计算两个正数的对数平均值。

    对数平均定义为 (b - a) / (ln(b) - ln(a))，适用于 a, b > 0。
    当 a = b 时，极限值为 a 本身。

    参数:
        a: 第一个正数
        b: 第二个正数

    返回:
        对数平均值

    异常:
        ValueError: 如果 a 或 b 不是正数
    """
    if a <= 0 or b <= 0:
        raise ValueError("对数平均要求输入为正数")

    if a == b:
        return a

    # 防止数值问题，确保计算稳定
    # 如果 a 和 b 非常接近，可以用泰勒展开，但这里直接计算
    return (b - a) / math.log(b / a)


def arlog_geo_mean(a: float, b: float) -> float:
    """
    计算算术平均与对数平均的几何平均，即 sqrt(算术平均 * 对数平均)。

    该组合平均结合了算术平均和对数平均的特性，常用于某些工程或数学语境。
    要求 a, b > 0（因为对数平均需要正数）。

    参数:
        a: 第一个正数
        b: 第二个正数

    返回:
        算术平均与对数平均的几何平均值

    异常:
        ValueError: 如果 a 或 b 不是正数
    """
    if a <= 0 or b <= 0:
        raise ValueError("arlog_geo_mean 要求输入为正数")
    arith = arithmetic_mean(a, b)
    log = logarithmic_mean(a, b)
    return (arith * log) ** 0.5


# 示例测试
if __name__ == "__main__":
    # 测试算术平均
    print(f"算术平均(3, 5) = {arithmetic_mean(3, 5)}")
    # 测试对数平均
    print(f"对数平均(3, 5) = {logarithmic_mean(3, 5)}")
    print(f"对数平均(5, 5) = {logarithmic_mean(5, 5)}")
    # 测试组合平均
    print(f"arlog_geo_mean(3, 5) = {arlog_geo_mean(3, 5)}")
