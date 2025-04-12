from matplotlib import pyplot as plt
from base_functions import *

def step_line(m_eq, Y_H, Y_K, X_H, X_K):
    '''
    "Графическое" нахождения числа тепор ступеней
    '''
    a_rab_line = (Y_H - Y_K) / (X_K - X_H)
    b_rab_line = Y_K - a_rab_line * X_H
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
    print(a_rab_line, b_rab_line)

    print(full_steps + step_part)

    plt.figure()
    X_rab_left = X_step_mas[-1] - abs(X_step_mas[-1]) * dec(0.3)
    plt.plot((X_rab_left, X_K), (a_rab_line * X_rab_left + b_rab_line, Y_H), color=(0, 0.6902, 0.3137))  # рабочая
    plt.plot((X_H, X_K), (Y_K, Y_H), 'g.')
    plt.plot((X_H, X_K), (m_eq * X_H, m_eq * X_K), color=(0, 0.4392, 0.7529))  # равновесная
    plt.plot((X_H, X_K), (m_eq * X_H, m_eq * X_K), 'b.')
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
    plt.plot((X_H, X_K), (Y_K, Y_H), 'g.')
    plt.plot((X_H, X_K), (m_eq * X_H, m_eq * X_K), color=(0, 0.4392, 0.7529))  # равновесная
    plt.plot((X_H, X_K), (m_eq * X_H, m_eq * X_K), 'b.')
    # ax.quiver(pos_x, pos_y, u / norm, v / norm, angles="xy", zorder=5, pivot="mid")
    plt.plot(X_step_mas, Y_step_mas, '-r.')
    plt.plot((X_H, X_H), (Y_K, m_eq * X_H), 'k')
    plt.grid(visible=True, which='major', color=(0.3, 0.3, 0.3), linestyle='-')  # сетка
    plt.grid(visible=True, which='minor', color=(0.63, 0.63, 0.63), linestyle='-', alpha=0.2)  # сетка
    plt.minorticks_on()  # сетка
    plt.xlabel("X", fontsize=14, x=1)
    plt.ylabel("Y", fontsize=14, rotation='horizontal', y=1)
    plt.title(f"Равновесная и рабочая линии")
    plt.show()


def conc_transfer(z_indef, z_type, M_1, M_2, is_liquid, ind=''):
    '''
    Перевод между различными концентрациями. Пока есть между относительными, абсолютными, мольными и массовыми
    '''
    ind = ind if ind == '' else '_' + ind

    conc_names = ('z', 'Z', 'Z_mass', 'z_mass')
    intermediate_dict = {key: 0 for key in conc_names}
    intermediate_dict[z_type] = z_indef
    z, Z, Z_mass, z_mass = (intermediate_dict[key] for key in conc_names)

    carrier_name, part_name = ('L', 'x') if is_liquid else ('G', 'y')

    def from_z():
        nonlocal Z
        Z = z / (1 - z)
        print(rawTOfinal(f'{part_name.upper()}{ind} = {part_name}{ind} / (1 - {part_name}{ind}) = {z} / (1 - {z}) = {Z}'))
    def from_Z():
        if action == -1:
            nonlocal z
            z = Z / (Z + 1)
            print(rawTOfinal(f'{part_name}{ind} = {part_name.upper()}{ind} / ({part_name.upper()}{ind} + 1) = {Z} / ({Z} + 1) = {z}'))
        else:
            nonlocal Z_mass
            Z_mass = Z * M_1 / M_2
            print(rawTOfinal(f'{part_name.upper()} _mass {ind} = {part_name.upper()}{ind} * @@ M_A / M_{carrier_name} = {Z} * @@ {M_1} / {M_2} = {Z_mass}'))
    def from_Z_mass():
        if action == -1:
            nonlocal Z
            Z = Z_mass * M_2 / M_1
            print(rawTOfinal(f'{part_name.upper()}{ind} = {part_name.upper()} _mass {ind} * @@ M_{carrier_name} / M_A = {Z_mass} * @@ {M_2} / {M_1} = {Z}'))
        else:
            nonlocal z_mass
            z_mass = Z_mass / (Z_mass + 1)
            print(rawTOfinal(f'{part_name} _mass {ind} = {part_name.upper()} _mass {ind} / ({part_name.upper()} _mass {ind} + 1) = {Z_mass} / ({Z_mass} + 1) = {z_mass}'))
    def from_z_mass():
        nonlocal Z_mass
        Z_mass = z_mass / (1 - z_mass)
        print(rawTOfinal(f'{part_name.upper()} _mass {ind} = {part_name} _mass {ind} / (1 - {part_name} _mass {ind}) = {z_mass} / (1 - {z_mass}) = {Z_mass}'))


    func_names = {'z': from_z,
                  'Z': from_Z,
                  'Z_mass': from_Z_mass,
                  'z_mass': from_z_mass}
    actions_names = {'z': (1,),
                     'Z': (-1, 1),
                     'Z_mass': (-1, 1),
                     'z_mass': (-1,)
                     }
    for action in actions_names[z_type]:
        func_names[z_type]()
        name_num = conc_names.index(z_type) + action
        while 1 <= name_num < len(conc_names) - 1:
            next_name = conc_names[name_num]
            func_names[next_name]()
            name_num += 1


def driving_force(m_eq, Y_H, Y_K, X_H, X_K, for_gas=True, for_mol=True):
    '''
    Движущая сила массопередачи, только для прямой равновесной линии
    Можно сделать как по жидкой, так и по газовой фазам
    Автоопределение порядка вычитания, из рабочих равновесных или наоборот
    '''
    if for_mol:
        overline = ''
    else:
        overline = '_mass'
    is_sorbtion = Y_H > m_eq * X_K

    if for_gas:
        Y_eq_X_H = m_eq * X_H
        Y_eq_X_K = m_eq * X_K
        Delta_Y_verh = Y_K - Y_eq_X_H if is_sorbtion else Y_eq_X_H - Y_K
        Delta_Y_niz = Y_H - Y_eq_X_K if is_sorbtion else Y_eq_X_K - Y_H

        Y_eq_X_H_form = f'Y {overline} ^*_X {overline} _H = m {overline} * X {overline} _H = {m_eq} * {X_H} = {Y_eq_X_H}'
        Y_eq_X_K_form = f'Y {overline} ^*_X {overline} _K = m {overline} * X {overline} _K = {m_eq} * {X_K} = {Y_eq_X_K}'

        Delta_Y_verh_form = f' Delta Y {overline} _верх = Y {overline} _K - Y {overline} ^*_X {overline} _H = {Y_K} - {Y_eq_X_H} = {Delta_Y_verh}' if is_sorbtion\
                       else f' Delta Y {overline} _верх = Y {overline} ^*_X {overline} _H - Y {overline} _K = {Y_eq_X_H} - {Y_K} = {Delta_Y_verh}'
        Delta_Y_niz_form = f' Delta Y {overline} _низ = Y {overline} _H - Y {overline} ^*_X {overline} _K = {Y_H} - {Y_eq_X_K} = {Delta_Y_niz}' if is_sorbtion\
                      else f' Delta Y {overline} _низ = Y {overline} ^*_X {overline} _K - Y {overline} _H = {Y_eq_X_K} - {Y_H} = {Delta_Y_niz}'

        Delta_Y_sr = (Delta_Y_niz - Delta_Y_verh) / ln(Delta_Y_niz / Delta_Y_verh)
        Delta_Y_sr_form = (f' Delta Y {overline} _ср = ( Delta Y {overline} _низ - Delta Y {overline} _верх) / ln( Delta Y {overline} _низ / Delta Y {overline} _верх)'
                           f' = ({Delta_Y_niz} - {Delta_Y_verh}) / ln({Delta_Y_niz} / {Delta_Y_verh}) = {Delta_Y_sr}')
        answ_form = '\n'.join((Y_eq_X_H_form, Y_eq_X_K_form, Delta_Y_verh_form, Delta_Y_niz_form, Delta_Y_sr_form))
        return Delta_Y_sr, answ_form

    else:
        X_eq_Y_H = Y_H / m_eq
        X_eq_Y_K = Y_K / m_eq
        Delta_X_verh = X_K - X_eq_Y_H if not is_sorbtion else X_eq_Y_H - X_K
        Delta_X_niz = X_H - X_eq_Y_K if not is_sorbtion else X_eq_Y_K - X_H

        X_eq_Y_H_form = f'X {overline} ^*_Y {overline} _H =  Y {overline} _H / m {overline} = {m_eq} * {Y_H} = {X_eq_Y_H}'
        X_eq_Y_K_form = f'X {overline} ^*_Y {overline} _K = Y {overline} _K / m {overline} = {m_eq} * {Y_K} = {X_eq_Y_K}'

        Delta_X_verh_form = f' Delta X {overline} _верх = X {overline} _K - X {overline} ^*_Y {overline} _H = {X_K} - {X_eq_Y_H} = {Delta_X_verh}' if not is_sorbtion \
                       else f' Delta X {overline} _верх = X {overline} ^*_Y {overline} _H - X {overline} _K = {X_eq_Y_H} - {X_K} = {Delta_X_verh}'
        Delta_X_niz_form = f' Delta X {overline} _низ = X {overline} _H - X {overline} ^*_Y {overline} _K = {X_H} - {X_eq_Y_K} = {Delta_X_niz}' if not is_sorbtion \
                      else f' Delta X {overline} _низ = X {overline} ^*_Y {overline} _K - X {overline} _H = {X_eq_Y_K} - {X_H} = {Delta_X_niz}'

        Delta_X_sr = (Delta_X_niz - Delta_X_verh) / ln(Delta_X_niz / Delta_X_verh)
        Delta_X_sr_form = (
            f' Delta X {overline} _ср = ( Delta X {overline} _низ - Delta X {overline} _верх) / ln( Delta X {overline} _низ / Delta X {overline} _верх)'
            f' = ({Delta_X_niz} - {Delta_X_verh}) / ln({Delta_X_niz} / {Delta_X_verh}) = {Delta_X_sr}')
        answ_form = '\n'.join((X_eq_Y_H_form, X_eq_Y_K_form, Delta_X_verh_form, Delta_X_niz_form, Delta_X_sr_form))
        return Delta_X_sr, answ_form

# m_eq, Y_H, Y_K, X_H, X_K = map(dec, (0.527, 0.05263, 0.007895, 0, 0.06657))
# step_line(m_eq, Y_H, Y_K, X_H, X_K)
# # conc_transfer(dec(0.05782), 'z', 34, 29, False, 'H')