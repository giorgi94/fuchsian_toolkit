from fuchsian_toolkit import get_fuchsian_transform


def example_1():
    A_list = [
        Matrix([[S("-2+1/4"), S("I")], [S("0"), S("-1+2/4")]]),
        Matrix([[S("-1/4"), S("1/2")], [S("0"), S("-1/3")]]),
        Matrix([[S("1+1/4"), S("-I")], [S("0"), S("1/3")]]),
        Matrix([[S("1-1/4"), S("-1/2")], [S("0"), S("1-2/4")]]),
    ]

    s_list = [
        S("(-1/2)+(3/4)*I"),
        S("(1/3)+(3/2)*I"),
        S("(1/2)-(1/4)*I"),
        S("(-1/3)-(1/3)*I"),
    ]

    mu_list = [
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
    ]

    T, fuch = get_fuchsian_transform(A_list, s_list, mu_list)

    for A in fuch.A_list:
        print(A)

    print()

    print(fuch.C)

    print()

    print(T)


def example_2():
    A_list = [
        Matrix([[S("2+1/2"), S("0")], [S("-1/4"), S("1+1/4")]]),
        Matrix([[S("-1/2"), S("0")], [S("0"), S("-1+1/3")]]),
        Matrix([[S("-1+1/2"), S("0")], [S("1/4"), S("-1-1/3")]]),
        Matrix([[S("-1-1/2"), S("0")], [S("0"), S("1-1/4")]]),
    ]

    s_list = [
        S("-1/2"),
        S("I/3"),
        S("1/4"),
        S("-(5/2)*I"),
    ]

    mu_list = [
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
    ]

    T, fuch = get_fuchsian_transform(A_list, s_list, mu_list)

    for A in fuch.A_list:
        print(A)

    print()

    print(fuch.C)

    print()

    print(T)


def example_3():
    A_list = [
        Matrix([[S("-4"), S("-1")], [S("2"), S("-7/6")]]),
        Matrix([[S("4"), S("-14/3")], [S("1/2"), S("11/12")]]),
        Matrix([[S("1/4"), S("1/6")], [S("-9/4"), S("-1")]]),
        Matrix([[S("-1/4"), S("11/2")], [S("-1/4"), S("5/4")]]),
    ]

    s_list = [
        S("1/4+1/4*I"),
        S("-1/3+2/5*I"),
        S("-1/4-1/2*I"),
        S("1/10-3*I"),
    ]

    mu_list = [
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
    ]

    T, fuch = get_fuchsian_transform(A_list, s_list, mu_list)

    for A in fuch.A_list:
        print(A)

    print()

    print(fuch.C)

    print()

    print(T)


if __name__ == "__main__":
    # example_0()
    # example_1()
    # example_2()
    example_3()
