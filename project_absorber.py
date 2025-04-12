from absorbtion_func import driving_force
from appr_module import *
from diffusion import diff_in_gas, diff_in_liq, diff_in_diluted_solv
from props_for_abs import *  # селексол только для примера
from column_function import *
from equil_data import *
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


# кафедральный пример
# выбор веществ, указать по номеру в словаре
sub_A = {1: 'сероводород', 2: 'углекислый газ'}[1]  # поглощаемое

sub_L = {1: 'пропиленкарбонат', 2: 'метанол',
         3: 'селексол', 4: 'вода',
         5: 'флотореагент т-66', 6: 'метилпирролидон'}[3]  # поглотитель

sub_G = {1: 'природный газ', 2: 'метан',
         3: 'водород'}[1]  # инерт

p_1 = dec('7') * 10 ** 6  # давление в абсорбере
p_2 = dec('0.1') * 10 ** 6  # давление в десорбере
t_1 = 25  # температура в абсорбере
t_2 = 75  # температура в десорбере
y_H = dec('0.09')  # мольная доля поглощаемого в инерте (равна объёмной доле)
is_normal: bool = True  # расход в нормальных кубометрах
V_yH = dec('75_000') / 3600  # объёмный расход
phi = dec('0.9')  # стпень поглощения
r = dec('1.4')  # коэффициент избытка поглотителя

element_name, order_type, sizes = 'керамические кольца Рашига', 'нерегулярная', '25_25_3'
# размеры через _

Psi = dec('0.7')  # коэффициент смачивания насадки
n_w = dec('0.8')  # отношение скорости газа к скорости захлёбывания




T_1, T_2 = dec(t_1 + 273), dec(t_2 + 273)

def get_equil_data(sub_A, sub_L, t, p):
    '''
    Функция подбора равновесных данных по заданным температуре и давлению. Адекватно работает только с табличными данными
    :param sub_A: название вещества поглощаемого вещества
    :param sub_L: название поглотителя
    :param t: рабочая температура
    :param p: рабочее давление
    :return: X_eq_nparray, Y_eq_nparray - равновесные данные
    '''

    def find_koef_a_p():
        '''
        Функция, подбирающая коэффициенты перевода размерностей для растворимости и давления из тех, что даны в
        справочнике, в СИ
        '''
        if alpha_dim == 'м³/м³':
            koef_alpha = M_L / (rho_L * V_0m)
            koef_alpha_formula = f"M_L / ( rho _L * V^0_m) = {M_L} / ({rho_L} * {V_0m})"
        elif alpha_dim == 'см³/г':
            koef_alpha = M_L / (1000 * V_0m)
            koef_alpha_formula = f"M_L / (1000 * V^0_m) = {M_L} / (1000 * {V_0m})"
        else:
            koef_alpha = dec(input(('неверная размерность концентрации, введи коэф. руками: ')))
            koef_alpha_formula = 'koef_alpha_formula'

        if p_dim == 'МПа':
            koef_p = 10 ** 6
        elif p_dim == 'кПа':
            koef_p = 10 ** 5
        elif p_dim == 'мм рт.ст.':
            koef_p = 133.32
        elif p_dim == 'ат' or p_dim == 'кг/см²':
            koef_p = 98100
        else:
            koef_p = dec(input(('неверная размерность давления, введи коэф. руками: ')))
        print(rawTOfinal(f"k_ alpha = {koef_alpha_formula} = {koef_alpha}"))
        return dec(koef_alpha), dec(koef_p), koef_alpha_formula

    t_tuple, alpha_p_tuples, alpha_dim, p_dim = mixtures_list[f'{sub_A}-{sub_L}']()
    koef_alpha, koef_p, koef_alpha_formula = find_koef_a_p()  # подбираем коэффициенты для конц. и давления

    def find_ref_data():
        '''
        Поиск справочных данных по растворимости газа
        '''
        # самый простой случай, когда данная температура совпала со справочной, просто формируем вывод
        if t in t_tuple:
            i = t_tuple.index(t)
            p_tuple, alpha_tuple = array_float_to_dec(alpha_p_tuples[i][1]), array_float_to_dec(alpha_p_tuples[i][0])
            columns = [[f'P({sub_A}),', f"{sym['alpha']} при t={t} °C,"],
                       [f'{p_dim}', f"{alpha_dim}"]]
            table_data = []
            for j in range(len(p_tuple)):
                table_data.append(create_table_line(p_tuple[j], alpha_tuple[j]))
            print_table(table_data, columns, 1, max_width=150)
            return p_tuple, alpha_tuple

        # нетабличная температура, ищем интервал
        try:
            t_x1, t_x2, alpha_p_y1_tuple, alpha_p_y2_tuple = interval_finder(t, t_tuple, alpha_p_tuples)
            # можно интерполировать по температуре, т.к. совпало число известных конц. для разных температур И одинаковые давления
            if len(alpha_p_y1_tuple[0]) == len(alpha_p_y2_tuple[0]) and alpha_p_y1_tuple[1] == alpha_p_y2_tuple[1]:
                p_tuple = array_float_to_dec(alpha_p_y1_tuple[1])
                alpha_tuple = []
                for i in range(len(alpha_p_y1_tuple[0])):
                    alpha_tuple.append(
                        interpol(t, (t_x1, t_x2), (alpha_p_y1_tuple[0][i], alpha_p_y2_tuple[0][i]), f'alpha_{i}'))
                alpha_tuple = array_float_to_dec(alpha_tuple)

                # собираем части таблицы после интерполяции, шапка с тремя разными температурами
                columns = [
                    [f"t={t_x1} {dimensions['t']}", '', f"t={t_x2} {dimensions['t']}", '', f"t={t} {dimensions['t']}",
                     ''],
                    [f"P, {p_dim}", f"{sym['alpha']}, {alpha_dim}", f"P, {p_dim}", f"{sym['alpha']}, {alpha_dim}",
                     f"P, {p_dim}",
                     f"{sym['alpha']}, {alpha_dim}"]]
                # собираем данные для таблицы, если температура не справочная
                table_data = []
                for j in range(len(p_tuple)):
                    table_data.append(create_table_line(alpha_p_y1_tuple[1][j], alpha_p_y1_tuple[0][j],
                                                        alpha_p_y2_tuple[1][j], alpha_p_y2_tuple[0][j],
                                                        p_tuple[j], alpha_tuple[j]))
                print()
                print_table(table_data, columns, 1, max_width=150)
                return p_tuple, alpha_tuple

            # температура не совпала и не проинтерполировать по температуре
            elif len(alpha_p_y1_tuple[0]) != len(alpha_p_y2_tuple[0]):
                print('массивы растворимостей не равны, нужно найти самый короткий')
            else:
                print('массивы давлений не равны, двухмерная интерполяция')

        except:
            print('интерполировать нечего, ручной ввод')
        # сюда попадаем только если не произошла автоматическая интерполяция
        alpha_tuple = input('Введи концентрации (альфа) через пробел:').split()
        p_tuple = input('Введи давления через пробел:').split()
        alpha_tuple = array_float_to_dec(alpha_tuple)
        p_tuple = array_float_to_dec(p_tuple)
        return p_tuple, alpha_tuple

    p_tuple, alpha_tuple = find_ref_data()
    # собираем таблицу равновесных данных в относительных долях
    X_eq_list, Y_eq_list = [], []
    p_tuple_inPa = tuple(p_i * koef_p for p_i in p_tuple)
    columns = [["Pобщ, МПа", "P, МПа", f"{sym['alpha']}, {alpha_dim}", "X", "Y"]]
    table_data = []
    for i in range(len(p_tuple)):
        p_try = p_tuple_inPa[i]

        # общее давление больше парциального давления вещества и ненулевая растворимость
        if p > p_try:
            if alpha_tuple[-1] != 0:
                X_eq_list.append(dec(alpha_tuple[i] * koef_alpha))
                Y_eq_list.append(p_tuple_inPa[i] / (p - p_tuple_inPa[i]))
            else:
                X_eq_list.append(0)
                Y_eq_list.append(0)

        # общее давление равно или меньше парциального вещества (Y не имеет смысла, X будет определен при равенстве давлений)
        elif p <= p_try:
            end_ind = i
            for j in range(end_ind, len(p_tuple)):
                Y_eq_list.append(None)
            if p == p_tuple_inPa[end_ind]:
                X_app = alpha_tuple[end_ind] * koef_alpha
            else:
                pol_koef, eq_func, *_ = pol_koef_search(np.array(array_dec_to_float(p_i / 10 ** 6 for p_i in p_tuple_inPa)),
                                                        np.array(array_dec_to_float(alpha_tuple)),
                                                        approx_func_tuple, xlim= (0, p_tuple_inPa[-1] * dec(1.05) / 10 ** 6),
                                                        vertical_lines=(p / 10 ** 6,), data_name='Растворимость',
                                                        axis_names=('p', f'{sym["alpha"]}'))
                X_app = dec(eq_func(float(p / 10 ** 6), *pol_koef) * float(koef_alpha))
            for j in range(end_ind, len(p_tuple) + 1):
                X_eq_list.append(dec(X_app))
            break
    else:
        end_ind = len(p_tuple)


    for i in range(len(p_tuple)):
        table_data.append(create_table_line(p / 10 ** 6, p_tuple_inPa[i] / 10 ** 6, alpha_tuple[i], X_eq_list[i], Y_eq_list[i]))
        print(rawTOfinal(f"""
X_{i + 1} = alpha _{i + 1} * k_ alpha = {alpha_tuple[i]} * {koef_alpha} = {dec(X_eq_list[i])}
Y_{i + 1} = p_{i + 1} / (P_общ - p_{i + 1}) = ({p_tuple_inPa[i]}) / ({p} - {p_tuple_inPa[i]}) = {Y_eq_list[i]}"""))  # печать промежуточнх расчётов мольных отношений
    print()
    print_table(table_data, columns, 1, max_width=150)
    return np.array(array_dec_to_float(X_eq_list)), np.array(array_dec_to_float(Y_eq_list[0:end_ind]))  # данные
    # переводятся во флоат, чтобы можно было подогнать кривую через методы scipy


approx_func_tuple = tuple(Appr_func0(func_name) for func_name in ('линейная', 'квадратичная', 'кубическая'))
M_A, nu_A, mu_A, rho_A = sub_props(sub_A, 'A', t_1, T_1)
M_L, nu_L, mu_L, rho_L = sub_props(sub_L, 'L', t_1, T_1)
M_G, nu_G, mu_G, rho_G = sub_props(sub_G, 'G', t_1, T_1)
X_eq_nparray, Y_eq_nparray = get_equil_data(sub_A, sub_L, t_1, p_1)
X_eq_nparray = X_eq_nparray[0:len(Y_eq_nparray)]
# матбаланс
if not is_normal:
    V_0yH = (V_yH * T_0 * p_1) / (T_1 * p_0)
    print(rawTOfinal(f"V^0_yH = (V_yH * T_0 * p_1) / (T_1 * p_0) = ({V_yH} * {T_0} * {p_1}) / ({T_1} * {p_0}) = {V_yH}"))
else:
    V_0yH = V_yH
n_yH = V_0yH / V_0m
n_AH = n_yH * y_H
m_AH = n_AH * M_A
n_G = n_yH * (1 - y_H)
m_G = n_G * M_G
Delta_n_A = n_AH * phi
Delta_m_A = Delta_n_A * M_A
m_yH = m_G + m_AH
n_yK = n_yH - Delta_n_A
m_yK = m_yH - Delta_m_A
Y_H = y_H / (1 - y_H)
Y_K = Y_H * (1 - phi)
print(rawTOfinal(f'''
n_yH = V_yH / V^0_m = {V_yH} / {V_0m} = {n_yH}
n_AH = n_yH * y_H = {n_yH} * {y_H} = {n_AH}
m_AH = n_AH * M_A = {n_AH} * {M_A} = {m_AH}
n_G = n_yH * (1 - y_H) = {n_yH} * (1 - {y_H}) = {n_G}
m_G = n_G * M_G = {n_G} * {M_G} = {m_G}
 Delta _n_A = n_AH * phi = {n_AH} * {phi} = {Delta_n_A}
 Delta _m_A = Delta n_A * M_A = {Delta_n_A} * {M_A} = {Delta_m_A}
m_yH = m_G + m_AH = {m_G} + {m_AH} = {m_yH}
n_yK = n_yH - Delta n_A = {n_yH} - {Delta_n_A} = {n_yK}
m_yK = m_yH - Delta m_A = {m_yH} - {Delta_m_A} = {m_yK}
Y_H = y_H / (1 - y_H) = {y_H} / (1 - {y_H}) = {Y_H}
Y_K = Y_H * (1 - phi ) = {Y_H} * (1 - {phi}) = {Y_K}
'''))

X_eq_nparray_des, _ = get_equil_data(sub_A, sub_L, t_2, p_2)  # находим равновесные данные в условиях десорбера
X_H = dec(X_eq_nparray_des[0])
(pol_koef, eq_func,
X_eq_nparray_plot, Y_eq_nparray_plot) = pol_koef_search(X_eq_nparray, Y_eq_nparray,
                                                        approx_func_tuple,
                                                        xlim=(0, X_eq_nparray[-1] * 1.05), ylim=(0, Y_eq_nparray[-1] * 1.05),
                                                        vertical_lines=(), horizontal_lines=(Y_H, Y_K),
                                                        data_name='Равновесные данные', line_name='Равновесная линия',
                                                        axis_names=('X', 'Y'), spec_title= 'Рабочая линия будет лежать между красными линиями\nи левее точек')
X_eq_Y_H = dec(fsolve(lambda x: eq_func(x, *pol_koef) - float(Y_H), 0)[0])  # находим равновесную концентрацию
# в поглотителе на выходе, решая уравнение
n_Lmin = Delta_n_A / (X_eq_Y_H - X_H)
n_L = n_Lmin * r
m_L = n_L * M_L
X_K = X_H + Delta_n_A / n_L
n_xH = n_L + n_L * X_H
m_xH = m_L + n_L * X_H * M_A
n_xK = n_xH + Delta_n_A
m_xK = m_xH + Delta_m_A
a_rab_line = n_L / n_G
b_rab_line = Y_K - a_rab_line * X_H

print(rawTOfinal(f"""
Из@графика@равновесной@линии:
X^*@@(Y_H) = {X_eq_Y_H}
X_H = alpha * k_ alpha = {X_H}
n_Lmin = Delta n_A / (X^*@@(Y_H) - X_H) = {Delta_n_A} / ({X_eq_Y_H} - {X_H}) = {n_Lmin}
n_L = n_Lmin * r = {n_Lmin} * {r} = {n_L}
X_K = X_H + Delta n_A / n_L = {X_H} + {Delta_n_A} / {n_L} = {X_K}
n_xH = n_L + n_L * X_H = {n_L} + {n_L} * {X_H} = {n_xH}
m_xH = m_L + n_L * X_H * M_A = {m_L} + {n_L} * {X_H} * {M_A} = {m_xH}
n_xK = n_xH + Delta n_A = {n_xH} + {Delta_n_A} = {n_xK}
m_xK = m_xH + Delta m_A = {m_xH} + {Delta_m_A} = {m_xK}

Ответы@для@схемы:
Y_K = {Y_K}
n_yK = {n_yK}
m_yK = {m_yK}

Delta n_A = {Delta_n_A}
Delta m_A = {Delta_m_A}

n_G = {n_G}
m_G = {m_G}

Y_H = {Y_H}
n_yH = {n_yH}
m_yH = {m_yH}

X_H = {X_H}
n_xH = {n_xH}
m_xH = {m_xH}

n_L = {n_L}
m_L = {m_L}

X_K = {X_K}
n_xK = {n_xK}
m_xK = {m_xK}

Y_H = {Y_H}
X_K = {X_K}
Y_K = {Y_K}
X_H = {X_H}
a = n_L / n_G = {n_L} / {n_G} = {a_rab_line}
b = Y_K - a * X_H = {Y_K} - {a_rab_line} * {X_H} = {b_rab_line}
"""))

X_eq_list_plot, Y_eq_list_plot = [0, *X_eq_nparray_plot], [0, *Y_eq_nparray_plot]
plt.figure()
ax = plt.gca()
ax.set_xlim([0, X_K + X_K * dec(0.05)])
ax.set_ylim([0, Y_H + Y_H * dec(0.05)])
plt.plot((X_H, X_K), (Y_K, Y_H), color=(0, 0.6902, 0.3137))  # рабочая
plt.plot((X_H, X_K), (Y_K, Y_H), 'go')
plt.plot(X_eq_list_plot, Y_eq_list_plot, color=(0, 0.4392, 0.7529))  # равновесная
plt.plot(X_eq_nparray, Y_eq_nparray, 'r.')
plt.grid(visible=True, which='major', color=(0.3, 0.3, 0.3), linestyle='-')  # сетка
plt.grid(visible=True, which='minor', color=(0.63, 0.63, 0.63), linestyle='-', alpha=0.2)  # сетка
plt.minorticks_on()  # сетка
plt.xlabel("X", fontsize=14, x=1)
plt.ylabel("Y", fontsize=14, rotation='horizontal', y=1)
plt.title(f"Равновесная и рабочая линии")
plt.show()

# движущая сила
if len(pol_koef) == 1:
    Delta_Y_sr, Delta_Y_sr_form = driving_force(dec(pol_koef[0]), Y_H, Y_K, X_H, X_K, for_gas=True, for_mol=True)
    print(rawTOfinal(f'{Delta_Y_sr_form}'))
else:
    Y_nparray_Noy = np.linspace(float(Y_K), float(Y_H), 21)
    a_rab_line, b_rab_line = float(a_rab_line), float(b_rab_line)
    X_nparray_Noy = np.array([fsolve(lambda x: a_rab_line * x + b_rab_line - i, 0)[0] for i in Y_nparray_Noy])
    Y_eq_nparray_Noy = np.array([eq_func(i, *pol_koef) for i in X_nparray_Noy])
    Y_Yeq = np.array([1 / (Y_nparray_Noy[i] - Y_eq_nparray_Noy[i]) for i in range(len(Y_nparray_Noy))])

    plt.figure()
    ax = plt.gca()
    ord_max, ord_min = max(Y_Yeq), min(Y_Yeq)
    ord_min = ord_min - (ord_max - ord_min) / 8
    ord_max = ord_max + ord_max / 10
    ax.set_ylim([ord_min, ord_max])
    plt.plot((Y_nparray_Noy[0], Y_nparray_Noy[0]), (ord_min, Y_Yeq[0]), 'k--', lw=0.7)  # начало
    plt.plot((Y_nparray_Noy[-1], Y_nparray_Noy[-1]), (ord_min, Y_Yeq[-1]), 'k--', lw=0.7)  # конец
    plt.plot(Y_nparray_Noy, Y_Yeq, color=(0.035, 0.235, 0.57), marker='.')  # чеп
    plt.fill_between(Y_nparray_Noy, Y_Yeq, y2=0, fc='c')  # штриховка
    plt.grid(visible=True, which='major', color=(0.3, 0.3, 0.3), linestyle='-')  # сетка
    plt.grid(visible=True, which='minor', color=(0.63, 0.63, 0.63), linestyle='-', alpha=0.2)  # сетка
    plt.minorticks_on()  # сетка
    plt.xlabel("y", fontsize=12, x=1)
    plt.ylabel(r"$\frac{1}{Y - Y^*}$", fontsize=12, rotation='horizontal', y=1)
    plt.title('Число единиц переноса')
    plt.show()

    s_i_perenos = []
    for i in range(len(X_nparray_Noy) - 1):
        s_i_perenos.append((Y_Yeq[i] + Y_Yeq[i+1]) / 2 * (Y_nparray_Noy[i + 1] - Y_nparray_Noy[i]))

    N_Oy = dec(sum(s_i_perenos))
    s_i_perenos.append(N_Oy)
    Delta_Y_sr = (Y_H - Y_K) / N_Oy
    columns = [["X", "Y", "Y*", "1/(Y-Y*)", "si"],]
    table_data = []
    for i in range(len(X_nparray_Noy)):
        table_data.append(create_table_line(X_nparray_Noy[i], Y_nparray_Noy[i],
                                            Y_eq_nparray_Noy[i], Y_Yeq[i],
                                            s_i_perenos[i]))
    print_table(table_data, columns, 1, max_width=150)
    print(rawTOfinal(f"""
число@единиц@переноса@по@газовой@фазе:@N_Оy = {N_Oy}
 Delta Y_ср = (Y_H - Y_K) / N_Oy = ({Y_H} - {Y_K}) / {N_Oy} = {Delta_Y_sr}
"""))

# по жидкой фазе
    # pol_koef, eq_func, X_eq_plot, Y_eq_plot = pol_koef_search(Y_eq_mas_lit, X_eq_mas_lit)
    # X_eq_mas = np.array([eq_func(i, *pol_koef) for i in Y_mas])
    # X_Xeq = np.array([1 / (X_eq_mas[i] - X_mas[i]) for i in range(len(X_mas))])
    # print(Y_mas, X_mas, Y_eq_mas)
    # plt.figure()
    # plt.plot(Y_mas, Y_Yeq)
    # plt.show()
    # s_i_perenos = []
    # for i in range(len(X_mas) - 1):
    #     s_i_perenos.append((X_Xeq[i] + X_Xeq[i + 1]) / 2 * (X_mas[i + 1] - X_mas[i]))
    # s_i_perenos.append(0)
    #
    # plt.figure()
    # plt.plot(X_mas, X_Xeq)
    # plt.show()
    #
    # columns = [["Y", "X", "X*", "1/(X*-X)", "si"], ]
    # table_data = []
    # for i in range(len(X_mas)):
    #     table_data.append(list(map(lambda x: rawTOfinal(str(dec(x))),
    #                                (Y_mas[i], X_mas[i], X_eq_mas[i], X_Xeq[i], s_i_perenos[i]))))
    # print_table(table_data, columns, 1, max_width=150)
    # print(f"число единиц переноса по жидкой фазе: N_ox = {sum(s_i_perenos)}")


# расчёт диаметра
M_yH = M_A * y_H + M_G * (1 - y_H)
rho_yH = (M_yH * T_0 * p_1) / (V_0m * T_1 * p_0)
# вытягиваем предельную скорость
A_nas, B_nas, a_nas, eps, d_e, l_pack_el = packages(element_name, order_type, sizes)
Theta = A_nas - B_nas * (m_xK / m_yH) ** dec(0.25) * (rho_yH / rho_L) ** dec(0.125)
Omega = (a_nas * rho_yH) / (g * eps ** 3 * rho_L) * (mu_L * 1000 / dec(1.0026)) ** dec(0.16)
w_ypr = (10 ** Theta / Omega) ** dec(1/2)
w_y = w_ypr * n_w
print(rawTOfinal(f"""
t_1 = {t_1}
 rho _L = {rho_L}
 mu _L = {mu_L * 1000}
M_y_H = M_A * y_H + M_G * (1 - y_H) = {M_A} * {y_H} + {M_G} * (1 - {y_H}) = {M_yH}
 rho _yH = (M_yH * T_0 * p_1) / (V^0_m * T_1 * p_0) = ({M_yH} * {T_0} * {p_1}) / ({V_0m} * {T_1} * {p_0}) = {rho_yH}
 Theta = A - B * (m_xK / m_yH) ** 0.25 * ( rho _yH / rho _L) ** 0.125 =
= {A_nas} - {B_nas} * ({m_xK} / {m_yH}) ** 0.25 * ({rho_yH} / {rho_L}) ** 0.125 = {Theta}
 Omega = (a * rho _yH) / (g * eps ** 3 * rho _L) * (( mu _L)/ 1.0026) ** 0.16 =
 = ({a_nas} * {rho_yH}) / ({g} * {eps} ** 3 * {rho_L}) * ( {mu_L * 1000} / 1.0026) ** 0.16 = {Omega}
w_ypr = √(10 ** {Theta} / {Omega}) = {w_ypr}
w_y = w_yпр * n = {w_ypr} * {n_w} = {w_y}
"""))
if is_normal:
    V_0yH = V_yH
    V_yH = (V_yH * T_1 * p_0) / (T_0 * p_1)
    print(rawTOfinal(f"V_yH = V^0_yH * (T / T_0) * (p_0 / p) = {V_0yH} * {T_1} / {T_0} * ({p_0}) / ({p_1}) = {V_yH}"))
S_or = V_yH / w_y
D_or = ((4 * S_or) / pi) ** dec(1/2)
print(rawTOfinal(f"""
S_ор = V_yH / w_y = {V_yH} / {w_y} = {S_or}
D_ор = √((4 * S_ор) / pi ) = √((4 * {S_or}) / {pi}) = {D_or}
"""))
D = column_diam(D_or)


# коэффициент массоотдачи в газе
y_K = Y_K / (1 + Y_K)
M_yK = M_A * y_K + M_G * (1 - y_K)
rho_yK = (M_yK * T_0 * p_1) / (V_0m * T_1 * p_0)
rho_ysr = (rho_yH + rho_yK) / 2
mu_yH = M_yH / ((M_A * y_H) / mu_A + (M_G * (1 - y_H) / mu_G))
mu_yK = M_yK / ((M_A * y_K) / mu_A + (M_G * (1 - y_K) / mu_G))
mu_ysr = (mu_yH + mu_yK) / 2
S = (pi * D ** 2) / 4
w_yH = V_yH / S
V_m = V_0m * (T_1 / T_0) * (p_0 / p_1)
V_yK = n_yK * V_m
w_yK = V_yK / S
w_ysr = (w_yH + w_yK) / 2
Re_y = (w_ysr * d_e * rho_ysr) / (eps * mu_ysr)
print(rawTOfinal(f"""
y_K = Y_K / (1 + Y_K) = {Y_K} / (1 + {Y_K}) = {y_K}
M_yK = M_A * y_K + M_G * (1 - y_K) = {M_A} * {y_K} + {M_G} * (1 - {y_K}) = {M_yK}
 rho _yK = (M_yK * T_0 * p_1) / (V^0_m * T_1 * p_0) = ({M_yK} * {T_0} * {p_1}) / ({V_0m} * {T_1} * {p_0}) = {rho_yK}
 rho _yср = ( rho _yH + rho _yK) / 2 = ({rho_yH} + {rho_yK}) / 2 = {rho_ysr}
 mu _yH = M_yH / ((M_A * y_H) / mu _A + (M_G * (1 - y_H) / mu _G)) = {M_yH} / (({M_A} * {y_H}) / {mu_A * 10 ** 6} + ({M_G} * (1 - {y_H}) / {mu_G * 10 ** 6})) = {mu_yH * 10 ** 6}
 mu _yK = M_yK / ((M_A * y_K) / mu _A + (M_G * (1 - y_K) / mu _G)) = {M_yK} / (({M_A} * {y_K}) / {mu_A * 10 ** 6} + ({M_G} * (1 - {y_K}) / {mu_G * 10 ** 6})) = {mu_yK * 10 ** 6}
 mu _yср = ( mu _yH + mu _yK) / 2 = ({mu_yH * 10 ** 6} + {mu_yK * 10 ** 6}) / 2 = {mu_ysr * 10 ** 6}
S = ( pi * D ** 2) / 4 = ({pi} * {D} ** 2) / 4 = {S}
w_yH = V_yH / S = {V_yH} / {S} = {w_yH}
V_m = V^0_m * T / T_0 * p_0 / p = {V_0m} * {T_1} / {T_0} * ({p_0}) / ({p_1}) = {V_m}
V_yK = n_yK * V_m = {n_yK} * {V_m} = {V_yK}
w_yK = V_yK / S = {V_yK} / {S} = {w_yK}
w_yср = (w_yH + w_yK) / 2 = ({w_yH} + {w_yK}) / 2 = {w_ysr}
Re_y = (w_yср * d_e * rho _yср) / ( eps * mu _yср) = ({w_ysr} * {d_e} * {rho_ysr}) / ({eps} * {mu_ysr}) = {Re_y}
"""))

D_y = diff_in_gas(M_A, M_G, t_1, p_1, nu_A, nu_G)
Pr_y = mu_ysr / (rho_ysr * D_y)
if l_pack_el == 0:  # неупорядоченная насадка
    a, b, c = map(dec, (0.407, 0.665, 0.33))
    Nu_y = a * Re_y ** b * Pr_y ** c
    nusselt_formula = f"Nu_y = {a} * Re_y ** {b} * Рr_y ** {c} = {a} * {Re_y} ** {b} * {Pr_y} ** {c} = {Nu_y}"
else:
    a, b, c = map(dec, (0.167, 0.74, 0.33))  # упорядоченная
    Nu_y = a * Re_y ** b * Pr_y ** c * (l_pack_el / d_e) ** dec(-0.47)
    nusselt_formula = f"Nu_y = {a} * Re_y ** {b} * Рr_y ** {c} * (l/d_э) ** -0.47 = {a} * {Re_y} ** {b} * {Pr_y} ** {c} * ({l_pack_el} / {d_e}) ** -0.47 = {Nu_y}"

beta_yV = (Nu_y * D_y) / d_e
beta_y = beta_yV / V_m
print(rawTOfinal(f"""
Рr_y = mu _yср / ( rho _yср * D_y) = {mu_ysr} / ({rho_ysr} * {D_y}) = {Pr_y}
{nusselt_formula}
 beta _yV = (Nu_y * D_y) / d_e = ({Nu_y} * {D_y}) / {d_e} = {beta_yV}
 beta _y = beta _yV / V_m = {beta_yV} / {V_m} = {beta_y}
"""))


# коэффициент массоотдачи в жидкости
U = m_L / (rho_L * S)
Re_x = (4 * U * rho_L) / (a_nas * mu_L)

*_, mu_L20, rho_L20 = sub_props(sub_L, 'L', 20, 293)
print(rawTOfinal(f"""
U = m_L / ( rho _L * S) = {m_L} / ({rho_L} * {S}) = {U}
Re_x = (4 * U * rho _L) / (a * mu _L) = (4 * {U} * {rho_L}) / ({a_nas} * {mu_L}) = {Re_x}
"""))

if t_1 < 0:
    D_x20, _ = diff_in_liq(M_A, M_L, t_1, nu_A, nu_L, rho_L, mu_L, mu_L20, rho_L20, sub_L, sub_A)
    print('Пересчёт для отрицательных температур')
    D_x = diff_in_diluted_solv(M_L, t_1, nu_A, mu_L, mu_L20, rho_L20, sub_L, sub_A)
else:
    D_x20, D_x = diff_in_liq(M_A, M_L, t_1, nu_A, nu_L, rho_L, mu_L, mu_L20, rho_L20, sub_L, sub_A)

Pr_x = mu_L / (rho_L * D_x)
a, b, c = map(dec, (0.0021, 0.75, 0.5))
Nu_x = a * Re_x ** b * Pr_x ** c
delta_pr = (mu_L ** 2 / (rho_L ** 2 * g)) ** dec(1/3)
beta_xV = (Nu_x * D_x) / delta_pr
beta_x = (beta_xV * rho_L) / M_L
print(rawTOfinal(f"""
Рr_x = mu _L / ( rho _L * D_x) = {mu_L} / ({rho_L} * {D_x}) = {Pr_x}
Nu_x = {a} * Re_x ** {b} * Рr_x ** {c} = {a} * ({Re_x}) ** {b} * {Pr_x} ** {c} = {Nu_x}
 delta _pr = ( mu _L ** 2 / ( rho _L ** 2 * g)) ** (1/3) = ({mu_L} ** 2 / ({rho_L} ** 2 * {g})) ** (1/3) = {delta_pr}
 beta _xV = (Nu_x * D_x) / delta _пр = ({Nu_x} * {D_x}) / {delta_pr} = {beta_xV}
 beta _x = ( beta _xV * rho _L) / M_L = ({beta_xV} * {rho_L}) / {M_L} = {beta_x}
"""))

# коэффициент распределения
# если равновесная прямая:
if len(pol_koef) == 1:
    m_eq = dec(pol_koef[0])
# если нет:
else:
    print('равновесная линия не является прямой, коэффициент распределения находим как среднее значение производной '
          'в нескольких точках')
    X_nparray_for_m = np.linspace(float(X_H), float(X_K), 9)[1:-1]
    pol_koef_for_der = [*pol_koef, 0]
    pol_for_der = np.poly1d(pol_koef_for_der)
    derivated_pol = np.polyder(pol_for_der)
    m_eq_mas = [derivated_pol(el) for el in X_nparray_for_m]
    m_eq = dec(np.mean(m_eq_mas))
    columns = [["X", "m"], ]
    table_data = []
    for i in range(len(X_nparray_for_m)):
        table_data.append(create_table_line(X_nparray_for_m[i], m_eq_mas[i]))
    print_table(table_data, columns, 1, max_width=150)

# коэффициент массопередачи и высота
K_y = (1 / beta_y + m_eq / beta_x) ** -1
F = Delta_n_A / (K_y * Delta_Y_sr)
H_nas = F / (a_nas * S * Psi)
print(rawTOfinal(f"""
K_y = (1 / beta _y + m / beta _x) ** -1 = (1 / {beta_y} + {m_eq} / {beta_x}) ** -1 = {K_y}
F = Delta n_A / (K_y * Delta Y_sr) = {Delta_n_A} / ({K_y} * {Delta_Y_sr}) = {F}
"""))
H_col = packed_column_sizes(D_or, F, a_nas, Psi)

print(rawTOfinal(f'''
V_yH = {V_yH * 3600 if not is_normal else V_0yH * 3600}
m_L = {m_L} = {m_L * dec(3.6)}
'''))

# ответы для таблицы
decimal.getcontext().prec = 4
for el in (Delta_n_A, n_G, Y_H, Y_K, X_H, m_eq, n_Lmin, n_L, X_K, Delta_Y_sr, w_ypr, D_or, int(Re_y), D_y * 10 ** 7, Pr_y, Nu_y, beta_y, Re_x, D_x20 * 10 ** 9, D_x * 10 ** 9, Pr_x, delta_pr * 10 ** 6, Nu_x, beta_x, K_y, F, H_nas, H_col):
    print(dTc(dec(el)))