from base_functions import *
p_0 = 101325
T_0 = 273


def diff_in_gas(M_A, M_G, t, p, nu_A, nu_G, sub_A=''):
    M_A, M_G, t, nu_A, nu_G = map(dec, (M_A, M_G, t, nu_A, nu_G))
    p = int(p)
    T = t + 273

    D_y_calc = (dec(4.22 * 10 ** (-2)) * T ** dec(1.5)) / (p * (nu_A ** dec(1 / 3) + nu_G ** dec(1 / 3)) ** 2) * (
                1 / M_A + 1 / M_G) ** dec(1 / 2)
    print(rawTOfinal(f'''
D_y^(расч) = (4.22 * 10 ** (-2) * T ** (1.5)) / (p * ( nu _A ** (1/3) + nu _G ** (1/3)) ** 2) * (1 / M_A + 1 / M_G) ** (1 / 2) = (4.22 * 10 ** (-2) * {T} ** (1.5)) / ({p} * ({nu_A} ** (1/3) + {nu_G} ** (1/3)) ** 2) * (1 / {M_A} + 1 / {M_G}) ** (1 / 2) = {D_y_calc} {dimensions['trans_k']}
'''))

    if sub_A != '':  # если диффундирующее вещество указано, то проводится расчёт по известным опытным данным
        # коэффициент диффузии вещества в воздухе при 273К, 0 означает отсутствие данных
        D_y_0 = {'азот': 13.2,
                 'оксид азота': 14.5,
                 'аммиак': 19.8,
                 'водород': 61.1,
                 'кислород': 17.8,
                 'сероводород': 12.7,
                 'диоксид серы': 12.2,
                 'углекислый газ': 13.8,
                 'угарнй газ': 20.2,
                 'хлор': 12.4,
                 'хлороводород': 13.0,
                 'ацетилен': 20.6,
                 'бутан': 0,
                 'метан': 22.3,
                 'пропан': 0,
                 'этан': 0,
                 'этилен': 15.2,
                 'пропилен': 0}
        D_y_0 = D_y_0.get(sub_A, 0)
        if D_y_0 != 0:
            D_y_0 = dec(D_y_0 * 10 ** (-6))
            D_y_exp = D_y_0 * dec(p_0 / p) * (T / T_0) ** dec(1.5)
            Delta = (D_y_exp - D_y_calc) / D_y_exp
            print(rawTOfinal(f"""
D_y^0 = {D_y_0} {dimensions['trans_k']}
D_y^(эксп) = D_y^0 * (p_0 / p) * (T / T_0) ** (1.5) = {D_y_0} * ({p_0} / {p}) * ({T} / {T_0}) ** (1.5) = {D_y_exp} {dimensions['trans_k']}
 Delta = (D_y^(эксп) - D_y^(расч)) / D_y^(эксп) = ({D_y_exp} - {D_y_calc}) / ({D_y_exp}) = {Delta} = {Delta * 100} %
"""))
        else:
           print('Данных по диффузии газа в воздухе нет, D_y^(эксп) не известен')

    return D_y_calc


def diff_in_liq(M_A, M_L, t, nu_A, nu_L, rho_L, mu_L, mu_L20, rho_L20, sub_L, sub_A, ind_for_D=''):
    A_diff_dict = {'этанол': dec(1.24), 'метанол': dec(1.19), 'уксусная кислота': dec(1.27)}
    B_diff_dict = {'этанол': dec(2.0), 'метанол': dec(2.0), 'ацетон': dec(1.15), 'вода': dec(4.7)}
    A_diff = A_diff_dict.get(sub_A, 1)
    B_diff = B_diff_dict.get(sub_L, 1)
    M_A, M_L, t, nu_A, nu_L, rho_L, mu_L = map(dec, (M_A, M_L, t, nu_A, nu_L, rho_L, mu_L ))
    D_x20 = dec(10 ** (-6)) / (A_diff * B_diff * (mu_L20 * 1000) ** dec(1 / 2) * (nu_A ** dec(1 / 3) + nu_L ** dec(1 / 3)) ** 2) * (1 / M_A + 1 / M_L) ** dec(1 / 2)
    b_temp_koef_diff = dec(0.2) * ((mu_L20 * 1000) ** dec(1 / 2)) / (rho_L20 ** dec(1 / 3))
    D_x = D_x20 * (1 + b_temp_koef_diff * (t - 20))
    print(rawTOfinal(f"""
D{ind_for_D}_x20 = 10 ** (-6) / (A * B * √( mu _L20) * ( nu _A ** (1/3) + nu _L ** (1/3)) ** 2) * √(1 / M_A + 1 / M_L) = 
= 10 ** (-6) / ({A_diff} * {B_diff} * √({mu_L20 * 1000}) * ({nu_A} ** (1/3) + {nu_L} ** (1/3)) ** 2) * √(1 / {M_A} + 1 / {M_L}) = {D_x20} {dimensions['trans_k']}
b = 0.2 * √( mu _L20) / ∛( rho _L20) = 0.2 * √({(mu_L20 * 1000)}) / ∛({rho_L20}) = {b_temp_koef_diff}
D{ind_for_D}_x = D{ind_for_D}_x20 * (1 + b * (t - 20)) = {D_x20} * (1 + {b_temp_koef_diff} * ({t} - 20)) = {D_x} {dimensions['trans_k']}
"""))

    return D_x20, D_x


def diff_in_water_exp(t, mu_L20, rho_L20, sub_A):
    # диффузия газов в воде, см. задачник Носырева 2 сем с. 124
    D_x_0 = {'азот': 1.9,
             'оксид азота': 1.8,
             'аммиак': 1.8,
             'водород': 5.3,
             'кислород': 2.1,
             'сероводород': 1.6,
             'углекислый газ': 1.8,
             'хлор': 1.6,
             'хлороводород': 2.3}  # при 12 по Цельсию
    D_x_0 = dec(D_x_0.get(sub_A, 0) * 10 ** (-9))
    if sub_A == 'хлороводород':
        print('для HCl дана диффузия при 12, нужно исправить программу')
    if D_x_0 == 0:
        print('экспериментальное значение диффузии газа в воде не известно')
    else:
        b_temp_koef_diff = dec(0.2) * ((mu_L20 * 1000) ** dec(1 / 2)) / (rho_L20 ** dec(1 / 3))
        D_x_exp = D_x_0 * (1 + b_temp_koef_diff * (t - 20))
        print(rawTOfinal(f"""
D_x^0 = {D_x_0} {dimensions['trans_k']}
D_x^(эксп) = D_x^0 * (1 + b * (t - 20)) = {D_x_0} * (1 + {b_temp_koef_diff} * ({t} - 20)) = {D_x_exp} {dimensions['trans_k']}
"""))
        return D_x_exp


def diff_in_water_calc_and_exp(M_A, M_L, t, nu_A, nu_L, rho_L, mu_L, mu_L20, rho_L20, sub_L, sub_A):

    _, D_x_calc = diff_in_liq(M_A, M_L, t, nu_A, nu_L, rho_L, mu_L, mu_L20, rho_L20, sub_L, sub_A, ind_for_D='^(расч) ')
    D_x_exp = diff_in_water_exp(t, mu_L20, rho_L20, sub_A)

    if D_x_calc is not None:
        Delta = (D_x_exp - D_x_calc) / (D_x_exp)
        print(rawTOfinal(f"""
 Delta = (D_x^(эксп) - D_x^(расч)) / D_x^(эксп) = ({D_x_exp} - {D_x_calc}) / ({D_x_exp}) = {Delta} = {Delta * 100} %
"""))


def diff_in_diluted_solv(M_L, t, nu_A, mu_L, mu_L20, rho_L20, sub_L, sub_A):
    beta = {'вода': 2.6,
            'метанол': 1.9,
            'этанол': 1.5}
    beta = dec(beta.get(sub_L, 1))

    k1, k2, k3 = dec(7.4), dec(10 ** (-12)), dec(0.6)
    D_x_calc = k1 * k2 * ((t + 273) * sqrt(beta * M_L)) / (mu_L * 1000 * nu_A ** k3)
    D_x_exp = diff_in_water_exp(t, mu_L20, rho_L20, sub_A)
    if D_x_calc is not None:
        Delta = (D_x_exp - D_x_calc) / (D_x_exp)
        print(rawTOfinal(f"""
D_x^(расч) = 7.4 * 10 ** (-12) * ((t + 273) * sqrt ( beta * M_L)) / ( mu _L * nu _A ** 0.6) = 7.4 * 10 ** (-12) * (({t} + 273) * sqrt ({beta} * {M_L})) / ({mu_L * 1000} * {nu_A} ** 0.6) = {D_x_calc}
 Delta = (D_x^(эксп) - D_x^(расч)) / D_x^(эксп) = ({D_x_exp} - {D_x_calc}) / ({D_x_exp}) = {Delta} = {Delta * 100} %
"""))
    return D_x_calc