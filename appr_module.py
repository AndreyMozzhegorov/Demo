import numpy as np
from scipy.optimize import curve_fit

from base_functions import *
from matplotlib import pyplot as plt

import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)

class Appr_func0:
    '''
    Аппроксимирующие функции. Экземпляр содержит три атрибута: название функции, саму функцию и формулу
    '''
    def __init__(self, func_name):
        self.name = func_name
        self.func = {'линейная': lambda x, a: a * x,
                     'квадратичная': lambda x, a, b: a * x ** 2 + b * x,
                     'кубическая': lambda x, a, b, c: a * x ** 3 + b * x ** 2 + c * x
                     }[func_name]
        self.formula = {'линейная': lambda args: f'Y^* = a * x = {args[0]} * X',
                        'квадратичная': lambda args: f'Y^* = a * X ** 2 + b * X = {args[0]} * X ** 2 + {args[1]} * X',
                        'кубическая': lambda args: f'Y^* = a * X ** 3 + b * X ** 2 + c * X = '
                                                      f'{args[0]} * X ** 3 + {args[1]} * X ** 2 + {args[2]} * X'
                        }[func_name]


def pol_koef_search(x_nparray, y_nparray, approx_func_tuple,
                    xlim=(), ylim=(),
                    vertical_lines=(), horizontal_lines=(),
                    data_name='Исходные данные', line_name='Линия аппроксимации',
                    axis_names=('x', 'y'), spec_title=''):
    '''
    Функция подбора уравнения аппроксимации
    :param x_nparray: массив равновесных X
    :param y_nparray: массив равновесных Y
    :return: pol_koef_tuple - коэф. подобранного уравнения, eq_func - вид подобранного уравнения,
     x_nparray_plot и y_nparray_plot - данные для построения равновесной линии (для графика с рабочей)
    '''

    print(f'Наносим {data_name.lower()} на график, чтобы была ясна зависимость между величинами - линейная или нет',
          'Красные линии обозначают интересующую область/точки', sep='\n')

    xlim = (x_nparray.min() * 0.9, x_nparray.max() * 1.05) if xlim == () else xlim
    ylim = (y_nparray.min() * 0.9, y_nparray.max() * 1.05) if ylim == () else ylim

    plt.figure()
    ax = plt.gca()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.plot(x_nparray, y_nparray, 'r.')  # Аппроксимируемые точки
    plt.plot((0, x_nparray[-1]), (0, y_nparray[-1]), 'k', lw=0.6)  # прямая для сравнения
    for line in vertical_lines:
        plt.plot((line, line), (ylim), 'r-', lw=0.6)  # вертикаль
    for line in horizontal_lines:
        plt.plot((xlim), (line, line), 'r-', lw=0.6)  # горизонталь
    plt.grid(visible=True, which='major', color=(0.3, 0.3, 0.3), linestyle='-')  # сетка
    plt.grid(visible=True, which='minor', color=(0.63, 0.63, 0.63), linestyle='-', alpha=0.2)  # сетка
    plt.minorticks_on()  # сетка
    plt.xlabel(axis_names[0], fontsize=14, x=1)
    plt.ylabel(axis_names[1], fontsize=14, rotation='horizontal', y=1)
    plt.title(f"{data_name}, выбор точек\n{spec_title}")
    plt.show()

    # выбор данных
    while True:
        recalc_interval = input('номера первой и последней точек, по которым аппроксимируем данные, через пробел два числа'
                                '\n(нумерация точек с единицы), 0, если берём все: ').split()
        if recalc_interval == ['0']:
            x_nparray_try, y_nparray_try = x_nparray, y_nparray
            break
        elif len(recalc_interval) == 2:
            try:
                start_p, fin_p = tuple(map(lambda x: int(x) - 1, recalc_interval))  # перевод номеров в индексы
                # уточняем интервал, если последняя точка после первой И первая не меньше первой из справочных И
                # последняя не выходит за пределы справочных
                if fin_p >= start_p and start_p >= 0 and fin_p <= len(x_nparray):
                    x_nparray_try, y_nparray_try = x_nparray[start_p:fin_p + 1], y_nparray[start_p:fin_p + 1]
                    if len(x_nparray_try) == 1:  # для случая одной точки оставляем возможность проложить прямую
                        x_nparray_try, y_nparray_try = np.array((0, x_nparray[start_p])), np.array((0, y_nparray[start_p]))
                    break
                else:
                    print('неверный ввод')
                    continue
            except:
                print('неверный ввод')
                continue
        print('неверный ввод')

    x_nparray_plot = np.linspace(x_nparray_try[0], x_nparray_try[-1], 100)
    for eq_func_class_exmpl in approx_func_tuple:
        try:
            print(f'Тип функции: {eq_func_class_exmpl.name}')
            pol_koef_tuple, *_ = curve_fit(eq_func_class_exmpl.func, x_nparray_try, y_nparray_try)
            y_nparray_plot = np.array([eq_func_class_exmpl.func(i, *pol_koef_tuple) for i in x_nparray_plot])

            plt.figure()
            ax = plt.gca()
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            plt.plot(x_nparray_plot, y_nparray_plot, color=(0, 0.4392, 0.7529))  # равновесная
            plt.plot(x_nparray_try, y_nparray_try, 'r.')
            plt.grid(visible=True, which='major', color=(0.3, 0.3, 0.3), linestyle='-')  # сетка
            plt.grid(visible=True, which='minor', color=(0.63, 0.63, 0.63), linestyle='-', alpha=0.2)  # сетка
            plt.minorticks_on()  # сетка
            plt.xlabel(axis_names[0], fontsize=14, x=1)
            plt.ylabel(axis_names[1], fontsize=14, rotation='horizontal', y=1)
            plt.title(f"{line_name}")
            plt.show()

            # подтверждение или дальнейший перебор
            end_flag = input('1 если аппроксимационная кривая подходит, 0 если перебираем дальше')
            # и защита от неверного ввода
            while end_flag not in ('1', '0'):
                print('неверный ввод')
                end_flag = input('1 если аппроксимационная кривая подходит, 0 если перебираем дальше')
            if end_flag == '1':
                print(rawTOfinal(eq_func_class_exmpl.formula(tuple(dec(koef) for koef in pol_koef_tuple)
                                                             )))
                return pol_koef_tuple, eq_func_class_exmpl.func, x_nparray_plot, y_nparray_plot
        except OptimizeWarning:
            print('мало точек для проведения аппроксимации\n')
            continue
    else:
        print('Испробованы все аппроксимирующие функции, повторяем снова\n')
        pol_koef_tuple, eq_func, x_nparray_plot, y_nparray_plot = pol_koef_search(x_nparray, y_nparray, approx_func_tuple, xlim, ylim, vertical_lines, horizontal_lines, data_name, line_name, axis_names)
    return pol_koef_tuple, eq_func, x_nparray_plot, y_nparray_plot