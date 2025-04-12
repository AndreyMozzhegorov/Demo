from base_functions import *
import matplotlib.pyplot as plt

# высота сегмента насадки
h_c = 3

def x_yeq_list_split(x_list, y_eq_list, x_W, x_F, x_P, y_eq_W, y_eq_F, y_eq_P):
    '''
    Функция разделения равновесных данных на верх и низ. Для ректификации
    '''
    x_rab_DOWN_list = (x_W,)
    y_eq_DOWN_list = (y_eq_W,)
    for i, el in enumerate(x_list):
        if el > x_W and el < x_F:
            x_rab_DOWN_list += (el,)
            y_eq_DOWN_list += (y_eq_list[i],)
    else:
        x_rab_DOWN_list += (x_F,)
        y_eq_DOWN_list += (y_eq_F,)

    x_rab_UP_list = (x_F,)
    y_eq_UP_list = (y_eq_DOWN_list[-1],)
    for i, el in enumerate(x_list):
        if el > x_F and el < x_P:
            x_rab_UP_list += (el,)
            y_eq_UP_list += (y_eq_list[i],)
    else:
        x_rab_UP_list += (x_P,)
        y_eq_UP_list += (y_eq_P,)
    x_dict_lists = {'full': x_list, 'down': x_rab_DOWN_list, 'up': x_rab_UP_list, 'rab': x_rab_DOWN_list + x_rab_UP_list[1:]}
    y_eq_dict_lists = {'full': y_eq_list, 'down': y_eq_DOWN_list, 'up': y_eq_UP_list, 'rab': y_eq_DOWN_list + y_eq_UP_list[1:]}
    return x_dict_lists, y_eq_dict_lists

def koef_for_liq_VEP(U_H, U_B, sizes, w_GW_by_w_pr):
    '''
    Коэффициенты для частной высоты единицы переноса по жидкости. Для ректификации
    '''
    c_dict = {'x': (31.62145, 35.04524, 38.37712, 42.12261, 45.8681, 49.38381, 52.07229, 55.10545, 58.29945, 60.75815, 63.14791, 65.9283, 68.40998, 70.70782, 72.77589, 75.09671, 77.48647, 80.74941, 82.17408, 83.69066, 85.04639, 86.6319, 87.43615, 88.14848),
              'y': (0.99199, 0.98381, 0.9782, 0.97095, 0.9595, 0.94174, 0.92724, 0.90855, 0.8819, 0.85993, 0.83492, 0.8057, 0.77321, 0.74165, 0.71103, 0.66569, 0.6194, 0.54835, 0.51188, 0.46817, 0.42423, 0.37397, 0.33868, 0.30876)}
    Phi_arrays = {'25_25_3': {'x': (1166.32335, 2304.95078, 4682.49952, 6700.88098, 10285.33523, 14057.78679, 18308.88969, 21852.00576, 29155.00261, 25605.79946, 35728.67318, 44802.11571, 55793.91146, 70446.88478, 86828.05139, 104828.23373, 121572.48052, 142292.96203, 170611.49338),
                              'y': (0.04123, 0.04177, 0.04332, 0.04508, 0.04786, 0.05076, 0.05364, 0.05622, 0.06075, 0.05886, 0.06596, 0.07203, 0.08007, 0.08994, 0.10258, 0.11659, 0.13082, 0.14747, 0.17444)},
                  '35_35_4': {'x': (1164.98431, 1802.59792, 2763.67214, 4276.27188, 6077.51616, 8268.54429, 13349.55634, 17709.0047, 24512.09669, 37927.93945, 48220.04178, 69084.50286, 103274.39946, 142948.29711, 174176.00787),
                              'y': (0.03389, 0.03701, 0.04109, 0.04561, 0.05023, 0.05485, 0.0636, 0.07003, 0.07922, 0.09438, 0.10564, 0.13036, 0.16663, 0.211, 0.25343)},
                  '50_50_5': {'x': (1126.81397, 1524.26578, 2196.3837, 2738.39056, 3691.53671, 5127.31706, 6677.82801, 9085.27681, 11900.84168, 15375.54408, 20608.56903, 29525.74684, 37624.17694, 53656.76636, 76873.65859, 119084.49796, 169439.60687),
                              'y': (0.05583, 0.05817, 0.06082, 0.06278, 0.06604, 0.06954, 0.07314, 0.07766, 0.08168, 0.08561, 0.09132, 0.09915, 0.10502, 0.11523, 0.12778, 0.15062, 0.17485)}
                  }
    c_for_VEP = interpol(w_GW_by_w_pr, c_dict['x'], c_dict['y'], need_to_print=False)
    Phi_H = interpol(U_H, Phi_arrays[sizes]['x'], Phi_arrays[sizes]['y'], need_to_print=False)
    Phi_B = interpol(U_B, Phi_arrays[sizes]['x'], Phi_arrays[sizes]['y'], need_to_print=False)
    print(f"""
Фактор орошения [1, с. 233, рис. 6.6] (см. Приложение, с. 86):
нижняя часть ФН = {Phi_H}, верхняя часть ФВ = {Phi_B}
Коэффициент скорости при {w_GW_by_w_pr} % и керамических кольцах размером
{sizes} мм [1, с. 233, рис. 6.6] (см. Приложение, с. 86): {c_for_VEP} 
""")
    return Phi_H, Phi_B, c_for_VEP

def koef_for_gas_VEP(w_GW_by_w_pr, sizes, D):
    '''
    Коэффициенты для частной высоты единицы переноса по жидкости. Для ректификации
    '''
    if D > dec(0.8):
        n_for_VEP = 1
        print('Показатель степени при диаметре более 800 мм: n = 1')
    else:
        n_for_VEP = dec(1.24)
        print('Показатель степени при диаметре до 800 мм: n = 1,24')

    psi_arrays = {'25_25_3': {'x': (0, 3.49581, 5.54089, 7.67789, 9.97574, 12.61826, 15.44461, 18.3399, 21.92454, 25.76194, 30.54146, 35.13716, 37.52692, 40.00859, 43.06473, 47.49957, 52.11824, 57.05862, 61.49346, 65.23895, 69.26018, 73.1895, 77.3486, 80.86431, 84.17321),
                              'y': (0, 17.91307, 30.34928, 41.72459, 52.68732, 62.29445, 72.8446, 80.9193, 88.994, 96.47931, 102.13749, 105.20234, 105.85067, 106.08643, 105.20234, 104.554, 101.90173, 98.18855, 93.82703, 89.40658, 84.8093, 79.79945, 74.7896, 69.72081, 64.71096)},
                  '35_35_4': {'x': (0, 2.82943, 3.9324, 4.94345, 5.79366, 6.82769, 8.62001, 9.81489, 11.51529, 13.14676, 15.26078, 17.48969, 19.97137, 23.28027, 27.20959, 31.39167, 35.13716, 39.22732, 43.59323, 48.09701, 52.37101, 57.24244, 62.27473, 65.42278, 75.57926),
                              'y': (0, 25.10368, 36.47898, 45.49671, 54.86808, 62.29445, 74.7896, 83.10006, 91.88203, 99.95673, 108.03143, 116.8134, 123.82719, 131.90189, 138.91568, 143.92553, 147.4619, 150.70357, 153.12008, 155.12402, 156.42069, 157.71736, 158.18887, 158.8372, 159.24978)},
                  '50_50_5': {'x': (0, 1.70349, 2.14008, 2.73752, 3.49581, 4.34601, 5.38004, 6.66684, 7.51704, 8.18342, 9.63106, 10.66509, 12.11273, 13.49144, 15.09993, 17.74246, 20.13222, 22.68283, 25.85386, 29.50743, 33.68951, 37.71075, 40.51412, 45.20173, 48.62552, 58.00073, 70),
                              'y': (0, 34.06247, 45.43777, 56.34156, 67.3043, 78.91536, 91.58733, 104.49507, 113.45386, 121.17492, 129.72113, 137.1475, 145.28114, 153.12008, 162.31463, 171.27342, 179.58388, 185.301, 191.60752, 196.85313, 201.68616, 204.75101, 206.04767, 207.81586, 208.69995, 208.69995, 208.69995)}
                         }
    psi_VEP = interpol(w_GW_by_w_pr, psi_arrays[sizes]['x'], psi_arrays[sizes]['y'], need_to_print=False)
    print(f'Коэффициент скорости при {w_GW_by_w_pr} % и керамических кольцах размером {sizes} мм [1, с. 233, рис. 6.6] (см. Приложение, с. 86): {psi_VEP}')
    return psi_VEP, n_for_VEP


def packages(element_name: str, order_type: str, sizes: str):
    """
    Функция для нахождения параметров насадки
    Возвращает A_nas, B_nas, a_nas, eps, d_e типа Decimal
    :param element_name: название элементов, поддерживает только "керамические кольца Рашига"
    :param order_type: тип укладки, "регулярная" или "нерегулярная"+
    :param sizes: размеры через подчёркивание в мм, дроби черкз точку, пример: "10_10_1.5"
    """

    def ceramic_Raschig():
        """
        Функция для нахождения параметров колец Рашига
        Возвращает A_nas, B_nas, a_nas, eps, d_e типа Decimal
        """
        def regular_packing(sizes):
            A_nas, B_nas = dec(0.062), dec(1.55)  # параметры насадки для предельной скорости, с. 197 Дытнерский, пакетная насадка
            size_50_50_5 = {'a_nas': 110, 'eps': 0.735, 'd_e': 0.027, 'l': 0.05}  # удельная площадь, свободный объём и экв. диаметр
            size_80_80_8 = {'a_nas': 80, 'eps': 0.72, 'd_e': 0.036, 'l': 0.08}
            size_100_100_10 = {'a_nas': 60, 'eps': 0.72, 'd_e': 0.048, 'l': 0.10}
            # словарь размеров
            sizes_dict = {'50_50_5': size_50_50_5,
                          '80_80_8': size_80_80_8,
                          '100_100_10': size_100_100_10}

            finded_parameters = sizes_dict[f'{sizes}']  # находим нужные параметры по словарю
            a_nas, eps, d_e, l_pack_el = map(dec, (finded_parameters['a_nas'], finded_parameters['eps'], finded_parameters['d_e'], finded_parameters['l']))
            return A_nas, B_nas, a_nas, eps, d_e, l_pack_el

        def irregular_packing(sizes):
            A_nas, B_nas = dec(-0.073), dec(1.75)  # параметры насадки для предельной скорости, с. 197 Дытнерский, внавал
            size_10_10_1dot5 = {'a_nas': 440, 'eps': 0.7, 'd_e': 0.006}  # удельная площадь, свободный объём и экв. диаметр
            size_15_15_2 = {'a_nas': 330, 'eps': 0.7, 'd_e': 0.009}
            size_25_25_3 = {'a_nas': 200, 'eps': 0.74, 'd_e': 0.015}
            size_35_35_4 = {'a_nas': 140, 'eps': 0.78, 'd_e': 0.022}
            size_50_50_5 = {'a_nas': 90, 'eps': 0.785, 'd_e': 0.035}
            # словарь размеров
            sizes_dict = {'10_10_1.5': size_10_10_1dot5,
                          '15_15_2': size_15_15_2,
                          '25_25_3': size_25_25_3,
                          '35_35_4': size_35_35_4,
                          '50_50_5': size_50_50_5}
            finded_parameters = sizes_dict[f'{sizes}']  # находим нужные параметры по словарю
            a_nas, eps, d_e = map(dec, (finded_parameters['a_nas'], finded_parameters['eps'], finded_parameters['d_e']))
            l_pack_el = 0
            return A_nas, B_nas, a_nas, eps, d_e, l_pack_el

        order_dict = {'регулярная': regular_packing,
                      'нерегулярная': irregular_packing}
        A_nas, B_nas, a_nas, eps, d_e, l_pack_el = order_dict[f'{order_type}'](sizes)  # находим параметры по типу укладки
        # и размерам элементов
        return A_nas, B_nas, a_nas, eps, d_e, l_pack_el

    packing_dict = {'керамические кольца Рашига': ceramic_Raschig}

    A_nas, B_nas, a_nas, eps, d_e, l_pack_el = packing_dict[f'{element_name}']()
    # если нужна формула для коэф. А, то принтуй выше при расчёте
    # вывод найденных значений насадки
    print(f"""Характеристики насадки (по Дытнерскому, с. 196):
    насадка {element_name} {sizes}
    коэффициент А для предельной скорости {A_nas}
    коэффициент В для предельной скорости {B_nas}
    удельная поверхность {a_nas}
    свободный объём (порозность) {eps}
    эквивалентный диаметр {d_e}
""")
    return A_nas, B_nas, a_nas, eps, d_e, l_pack_el

def height_of_parts(D):
    '''
    Поиск высот частей колонны, зависящих от её диаметра
    :return: z_B, z_H, h_p, table_of_heights - высоты верхней, нижней частей, распределительной тарелки и таблица этих высот
    '''
    if dec(0.4) <= D <= dec(1):
        z_B, z_H = dec('0.6'), dec('1.5')
    elif dec(1.2) <= D <= dec(2.2):
        z_B, z_H = dec('1'), dec('2')
    elif dec(2.4) <= D:
        z_B, z_H = dec('1.4'), dec('2.5')
    else:
        print("Ошибка в диаметре")
        return
    # ищем высоту распределительной секции по ключу (диаметру)
    if D > 2.8:
        h_p = dec(0.915)
    else:
        h_p = dec({'0.4': 0.185,
               '0.5': 0.215,
               '0.6': 0.315,
               '0.8': 0.350,
               '1': 0.470,
               '1.2': 0.510,
               '1.4': 0.520,
               '1.6': 0.645,
               '1.8': 0.705,
               '2': 0.730,
               '2.2': 0.745,
               '2.4': 0.845,
               '2.6': 0.900,
               '2.8': 0.915}[str(D)[0:3]])  # отбрасываем два нуля с конца чтоб совпала строка
    "высота перераспределительной тарелки"
    columns = [["D, м", "Z_в, м", "Z_н, м", "h_п, м"], ]
    table_data = [create_table_line(D, z_B, z_H, h_p), ]
    table_of_heights = (table_data, columns)
    return z_B, z_H, h_p, table_of_heights

def column_diam(D_or):
    '''
    Определение стандартного диаметра по ориентировочному
    '''
    D_mas = (0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.4, 7.0, 8.0, 9.0)
    for el in D_mas:
        if dec(el) >= D_or:
            print(f"стандартный диаметр колонны: {el} м")
            return  dec(el)
    else:
        print('Слишком большой ориентировочный диаметр, нет подходящих стандартных!')


def packed_column_sizes(D_or, F, a_nas, Psi, area='F'):
    '''
    Нахождение высоты и диаметра насадочной колонны
    :param D_or: Ориентировочный диаметр
    :param F: площадь массопередачи
    :param a_nas: удельная поверхность насадки
    :param Psi: коэффициент самчивания
    :param area: символ для площади в уравнении: F или A
    :return: H - высота колонны
    '''
    D = column_diam(D_or)
    S = (pi * D ** 2) / 4
    H_nas = F / (a_nas * S * Psi)
    N_c = math.ceil(H_nas / h_c)
    z_B, z_H, h_p, table_of_heights = height_of_parts(D)
    print_table(*table_of_heights, 1, max_width=150)
    H = z_H + N_c * h_c + (N_c - 1) * h_p + z_B
    print(rawTOfinal(f"""
S = ( pi * D ** 2) / 4 = ( {pi} * {D} ** 2) / 4 = {S} @{dimensions['S']}
H_нас = {area} / (a * S * Psi ) = {F} / ({a_nas} * {S} * {Psi}) = {H_nas} @{dimensions['l']}
N_c = H_нас / h_c = {H_nas} / {h_c} = {N_c}
H = z_H + N_c * h_c + (N_c - 1) * h_п + z_B = {z_H} + {N_c} * {h_c} + ({N_c} - 1) * {h_p} + {z_B} = {H} @{dimensions['l']}
"""))
    return H

def plate_column_size(D_or, h_PT, eta, m_eq, Y_H, Y_K, X_H, X_K, a_rab_line, b_rab_line):
    '''
    Нахождение высоты и диаметра тарельчатой колонны через ступенчатую линию
    :param D_or: Ориентировочный диаметр
    :param h_PT: Высота тарелки
    :param eta: КПД
    :param m_eq: коэффициент распределения
    :return: H - высота колонны
    '''
    D = column_diam(D_or)
    z_B, z_H, h_p, table_of_heights = height_of_parts(D)
    print_table(*table_of_heights, 1, max_width=150)

    Y_step_mas = [Y_H, m_eq * X_K]
    X_step_mas = [X_K, X_K]
    while Y_step_mas[-1] > Y_K:
        X_step_mas.append((Y_step_mas[-1] - b_rab_line) / a_rab_line)
        X_step_mas.append(X_step_mas[-1])
        Y_step_mas.append(Y_step_mas[-1])
        Y_step_mas.append(m_eq * X_step_mas[-1])
    Y_step_mas.append(Y_step_mas[-1])
    X_step_mas.append((Y_step_mas[-1] - b_rab_line) / a_rab_line)
    step_part = (X_H - X_step_mas[-2]) / (X_step_mas[-1] - X_step_mas[-2])
    full_steps = dec((len(X_step_mas) - 3) / 2)

    X_H, X_K, Y_K, Y_H, m_eq = map(float, (X_H, X_K, Y_K, Y_H, m_eq))

    # заменить функцие ступенчатой линии
    plt.figure()
    X_rab_left = X_step_mas[-1] - abs(X_step_mas[-1]) * dec(0.3)
    plt.plot((X_rab_left, X_K), (a_rab_line * X_rab_left + b_rab_line, Y_H), color=(0, 0.6902, 0.3137))  # рабочая
    plt.plot((X_H, X_K), (Y_K, Y_H), 'go')
    plt.plot((X_H, X_K), (m_eq * X_H, m_eq * X_K), color=(0, 0.4392, 0.7529))  # равновесная
    plt.plot((X_H, X_K), (m_eq * X_H, m_eq * X_K), 'bo')
    plt.grid(visible=True, which='major', color=(0.3, 0.3, 0.3), linestyle='-')  # сетка
    plt.grid(visible=True, which='minor', color=(0.63, 0.63, 0.63), linestyle='-', alpha=0.2)  # сетка
    plt.minorticks_on()  # сетка
    plt.xlabel("X", fontsize=14, x=1)
    plt.ylabel("Y", fontsize=14, rotation='horizontal', y=1)
    plt.title(f"Равновесная и рабочая линии")
    plt.show()

    X_rab_left = X_step_mas[-1] - abs(X_step_mas[-1]) * dec(0.3)
    plt.figure()
    plt.plot((X_rab_left, X_K), (a_rab_line * X_rab_left + b_rab_line, Y_H), color=(0, 0.6902, 0.3137))  # рабочая
    plt.plot((X_H, X_K), (Y_K, Y_H), 'go')
    plt.plot((X_H, X_K), (m_eq * X_H, m_eq * X_K), color=(0, 0.4392, 0.7529))  # равновесная
    plt.plot((X_H, X_K), (m_eq * X_H, m_eq * X_K), 'bo')
    # ax.quiver(pos_x, pos_y, u / norm, v / norm, angles="xy", zorder=5, pivot="mid")
    plt.plot(X_step_mas, Y_step_mas, '-ro')
    plt.plot((X_H, X_H), (Y_K, m_eq * X_H), 'k')
    plt.grid(visible=True, which='major', color=(0.3, 0.3, 0.3), linestyle='-')  # сетка
    plt.grid(visible=True, which='minor', color=(0.63, 0.63, 0.63), linestyle='-', alpha=0.2)  # сетка
    plt.minorticks_on()  # сетка
    plt.xlabel("X", fontsize=14, x=1)
    plt.ylabel("Y", fontsize=14, rotation='horizontal', y=1)
    plt.title(f"Равновесная и рабочая линии")
    plt.show()

    # plt.plot(x_mas_plot, y_mas_plot)
    # plt.plot(x_mas_plot, y_eq_mas_plot)
    # plt.plot(X_step_mas, Y_step_mas, '-ro')
    # plt.plot((X_H, X_K), (Y_K, Y_H), 'g*')
    # plt.grid(True, alpha=0.5)
    # plt.xlabel("X", fontsize=14)
    # plt.ylabel("Y", fontsize=14)
    # plt.show()

    N_TT = step_part + full_steps
    N_PT = math.ceil(N_TT / eta)
    H = z_H + (N_PT - 1) * h_PT + z_B
    print(rawTOfinal(f"""
N_TT = {step_part} + {full_steps} = {N_TT}
N_PT = N_TT / eta = {N_TT} / {eta} = {N_PT}
H = z_H + (N_PT - 1) * h_PT + z_B = {z_H} + ({N_PT} - 1) * {h_PT} + {z_B} = {H} @{dimensions['l']}
"""))
    pass
