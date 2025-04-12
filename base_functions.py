import decimal as decimal
import functools
import math as math


decimal.getcontext().prec = 4
decimal.getcontext().rounding = decimal.ROUND_HALF_UP

dimensions = {'rho': '@кг/м^3',
              'mu': '@мПа⋅с',
              'c': '@@Дж/(кг⋅К)',
              'la': '@Вт/(м⋅К)',
              'beta': '@К^(-1)',
              'sigma': '@мДж/м^2',
              'r': '@кДж/кг',
              'p': '@@Па',
              't': '@°C',
              '': '',
              'mu_g': '@мкПа⋅с',
              'Pr': '',
              'l': '@м',
              'S': '@м^2',
              'v': '@м/с',
              'V': '@м^3/с',
              'n_ps': '@кмоль/с',
              'm_ps': 'кг/с',
              'K': '@@Вт/(м^2⋅К)',
              'R_to': '@(м^2⋅К)/Вт',
              'trans_k': '@@м^2/с',
              'K_y': 'кмоль/(м^2⋅с)',
              'MM': 'кг/кмоль',
              'K_y_m': 'кг/(м^2⋅с)'
              }

sym = {'alpha': 'α',
       'beta': 'β',
       'gamma': 'γ',
       'delta': 'δ',
       'eps': 'ε',
       'dzeta': 'ζ',
       'eta': 'η',
       'theta': 'θ',
       'kappa': 'κ',
       'la': 'λ',
       'mu': 'μ',
       'nu': 'ν',
       'sigma': 'σ',
       'tau': 'τ',
       'phi': 'φ',
       'rho': 'ρ',
       'pi': 'π',
       'Delta': 'Δ',
       'Sigma': 'Σ',
       'Psi': 'Ψ',
       'Theta': 'Θ',
       'Omega': 'Ω',
       '**': '^',
       '*': '⋅',
       'x_mass_W': 'x̅_W',
       'x_mass_F': 'x̅_F',
       'x_mass_P': 'x̅_P',
       'x_mass_D': 'x̅_D',
       'x_mass_os': 'x̅_ос',
       'x_mass_susp': 'x̅_сусп',
       'x_mass_f': 'x̅_ф',
       'x_mass_pr': 'x̅_пр',
       'x_mass_': 'x̅_',
       'x_mass': 'x̅',
       '_mass': '̅',
       'MB.': '',
       'sub_BK.M': 'M_BK',
       'sub_HK.M': 'M_HK',
       'sqrt': '√',
       'tqrt': '∛',
       'qqrt': '∜',
       '_starred': '^*',
       'Pr': 'Рr'}  # Тут разные буквы Р



def dec(b, prec=4):
    '''
    Перевод числа в тип decimal
    :param b: преобразуемое число
    :param prec: число значимых цифр, по умолчанию 4
    '''
    if prec != 4:
        old_prec = decimal.getcontext().prec
        decimal.getcontext().prec = prec
    answ = decimal.Decimal(b) + 0
    if prec != 4:
        decimal.getcontext().prec = old_prec
    return answ


def prec5(func):
    '''
    Декоратор для установки точности расчёта функции до 5 значащих цифр
    '''
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        old_prec = decimal.getcontext().prec
        decimal.getcontext().prec = 5
        answ = func(*args, **kwargs)
        decimal.getcontext().prec = old_prec
        return answ
    return wrapper


def arith_mean(values:tuple, names:tuple, dim=''):
    '''
    Среднеарифметическое + формула
    порядок имён соответствует порядку значений, последнее имя -- название усреднённой перемнной
    '''

    mean_answ = sum((dec(val) / 2 for val in values))
    if '' not in names:
        str_for_join_sym = ' = ( ' + ' + '.join((name for name in names[:-1])) + f') / {len(values)} = '
    else:
        str_for_join_sym = ' = '
    str_for_join_num = '( ' + ' + '.join((str(value) for value in values)) + f') / {len(values)} = {mean_answ}{dimensions[dim]}'
    str_for_print = ' ' + names[-1] + str_for_join_sym + str_for_join_num
    print(rawTOfinal(str_for_print))
    return mean_answ

def array_dec_to_float(x_array):
    '''
    перевод чисел типа decimal из коллекции в кортеж чисел типа float
    '''
    answ = []
    for el in x_array:
        try:
            answ.append(float(el))
        except:
            answ.append(el)
    return tuple(answ)


def array_float_to_dec(x_array):
    '''
    перевод чисел типа float из коллекции в кортеж чисел типа decimal
    '''
    answ = []
    for el in x_array:
        try:
            answ.append(dec(el))
        except:
            answ.append(el)
    return tuple(answ)


def sign(value):
    '''
    Функция знака.
    :return: value < 0 -> -1; value >= 0 -> 1
    '''
    if value == 0:
        return 1
    return (None, 1, -1)[int(abs(value) / value)]


def str_sign(value):
    '''
    Знак для формулы-строки
    '''
    return (None, '+', '')[sign(value)]


# константы
pi = dec('3.142')
g = dec('9.81')
p_0 = 101325
T_0 = dec('273')
V_0m = dec('22.414')
e_Napier = dec('2.718')

# мат. функции для децимал
def lg(num):
    return dec(math.log10(num))
def ln(num):
    return dec(math.log(num))
def sqrt(num):
    return dec(math.sqrt(num))

def dTc(dot_form):
    ''' замена точек на запятые'''
    dot_form = str(dot_form)
    comma_form = ''.join([i if i != '.' else ',' for i in dot_form])
    return comma_form

def rawTOfinal(raw_form):
    '''
    подготовка формул-строк к печати
    '''
    raw_list = raw_form.split(' ')
    final = ''
    for el in raw_list:
        if el in sym:
            final += sym[el]
            continue
        answ = ''
        flag_exp = False
        for char in el:
            if flag_exp == True:
                flag_exp = False
                if char == '+':
                    continue
                else:
                    answ += char
            elif char == '@':
                answ += ' '
            elif char == ',':
                answ += '.'
            elif char == '.':
                answ += ','
            elif (char == 'E') and ('+' in el or '-' in el):
                answ += '⋅10^'
                flag_exp = True
            else:
                answ += char
        if ',' in answ and (answ[-2:] == '00' or answ[-3:] == '000') and answ[-3] != ',':
            answ = answ[:-2]
        final += answ
    return final

def interval_finder(x, x_mas, y_mas):
    '''
    Нахождение интервала (простой перебор)
    '''
    for i in range(len(x_mas)):
        # нормальный вывод
        if x_mas[i] <= x <= x_mas[i+1] or x_mas[i] >= x >= x_mas[i+1]:
            answ = (x_mas[i], x_mas[i + 1],
                    y_mas[i], y_mas[i + 1])
            return answ
    # вывод, если ничего не нашли
    else:
        print('не попали в интервал')
        answ = (x_mas[-2], x_mas[-1],
                y_mas[-2], y_mas[-1])
        return answ

def interpol(x, x_mas, y_mas, name='x', need_to_print=True, need_float = False):
    if not need_float:
        x1, x2, y1, y2 = map(dec, interval_finder(x, x_mas, y_mas))
        x = dec(x)
        y = (y1 + (y2 - y1) * (x - x1) / (x2 - x1))
        str_form = f'{name}={y1}+({y2}-{y1})⋅({x}-{x1})/({x2}-{x1})={y} '
        if need_to_print:
            print(rawTOfinal(str_form))
        return y
    else:
        x1, x2, y1, y2 = interval_finder(x, x_mas, y_mas)
        x, x1, x2, y1, y2 = array_dec_to_float((x1, x2, y1, y2))
        y = (y1 + (y2 - y1) * (x - x1) / (x2 - x1))
        str_form = f'{name}={y1}+({y2}-{y1})⋅({x}-{x1})/({x2}-{x1})={y} '
        if need_to_print:
            print(rawTOfinal(str_form))
        return y



def create_table_line(*args):
    """переводит кортеж данных в децимал, а затем в строку, отсекая лишние знаки
    каждый отдельный элемент отправляется затем в отдельную колонку
    возвращает список, который нужно поместить в список"""
    table_line = []
    for arg in args:
        try:
            table_line.append(dTc(str(dec(arg))))
        except:
            table_line.append(arg)
    return table_line


def print_table(data, columns, indent, max_width=100):
    # стырено и модифицировано https://tproger.ru/articles/kak-napechatat-tablicu-s-pomoshhju-f-string
    # data — список списков, данные таблицы
    # columns — список списков, названия колонок таблицы
    # indent — отступ от края колонки
    # max_widt – допустимая ширина таблицы
    # max_columns — список максимальной длинны строки колонок
    # max_columns_title — список максимальной ширины колонок шапки
    # width — список ширины каждой колонки таблицы для печати

    # расчёт макимальной ширины колонок таблицы
    max_columns = []
    for col in zip(*data):
        len_el = []
        [len_el.append(len(str(el))) for el in col]
        max_columns.append(max(len_el))

    # вычислить максимальную длинну колонки шапки таблицы
    max_columns_title = []
    for col in zip(*columns):
        max_columns_title.append(max([len(el) for el in col]))

    # печать таблицы
    for col in columns:
        width = []
        for n, c in enumerate(col):

            # сравниваем максимальную колонку шапки с макс колонкой таблицы
            if max_columns[n] >= max_columns_title[n]:
                w = max_columns[n] + indent
                width.append(w)
            else:
                w = max_columns_title[n] + indent
                width.append(w)

            # пишем название колонок в две строки
            if sum(width) <= max_width:
                print(f'{c}\t', end='')  # выравнивание по центру print(f'{c:^{w}}', end='')
            else:
                print('Ширина таблицы больше допустимого значения')
                return
        print()

    # печать разделителя шапки '='
    print(f"{'=' * (sum(width))}")

    # печать тела таблицы
    for el in data:
        for n, col in enumerate(el):
            print(f'{col}\t', end='')  # выравнвание по правому краю print(f'{col:>{width[n]}}', end='')
        print()



mol_mass_dict = {'вода': 18.015,
                 'водород': 2.016,
                 'азот': 28.01,
                 'аммиак': 17.03,
                 'диоксид серы': 64.066,
                 'углекислый газ': 44.01,
                 'сероводород': 34.082,
                 'анилин': 93.128,
                 'ацетон': 58.08,
                 'бензол': 78.114,
                 'бромбензол': 157.01,
                 'бутанол': 74.123,
                 'воздух': 28.97,
                 'гексан': 86.177,
                 'гептан': 100.204,
                 'дихлорэтан': 98.959,
                 'изопропанол': 60.096,
                 'о-ксилол': 106.167,
                 'п-ксилол': 106.167,
                 'метанол': 32.042,
                 'метилацетат': 74.079,
                 'муравьиная кислота': 46.026,
                 'нитробензол': 123.111,
                 'октан': 114.231,
                 'пропанол': 60.096,
                 'сероуглерод': 76.143,
                 'тетрахлорметан': 153.822,
                 'толуол': 92.141,
                 'уксусная кислота': 60.053,
                 'хлорбензол': 112.558,
                 'хлороформ': 119.377,
                 'циклогексан': 84.161,
                 'этанол': 46.069,
                 'этилацетат': 88.106,
                 'этилбензол': 106.167,
                 'диэтиловый эфир': 74.123,
                 'фенол': 94.11,
                 'этиленгликоль': 62.07,
                 'изопропилбензол': 120.2
                 }


def pressure_transfer(p_old, measure, to_Pa=True):
    '''
    Конвертер давления в Паскали либо из них
    :param p_old: старое давление
    :param measure: at, mm Hg, bar, MPa - размерность из или в которую переводить
    :param to_Pa: в Паскали или нет
    :return:
    '''
    koef = dec({'at': 98100,
            'mm Hg': 133.32,
            'bar': 10 ** 5,
            'MPa': 10 ** 6}[measure])
    if to_Pa:
        p_new = str(int(p_old * koef))
        old_prec = decimal.getcontext().prec
        decimal.getcontext().prec = len(p_new)
        p_new = dec(p_new)
        decimal.getcontext().prec = old_prec
        print(rawTOfinal(f'p = {p_new}@Па'))
    else:
        p_new = dec(p_old) / koef
        print(rawTOfinal(f'p = {p_new}@{measure}'))
    return p_new