pkg load symbolic
syms t Vcc tau R_time_delay C_time_delay

% RC-цепочка на CONT
R_time_delay_v = 2.4e3; C_time_delay_v = 470e-6 + 470e-6;
Vcc_v = 12;

tau_time_delay = R_time_delay * C_time_delay
time_delay_factor = 5

% Построение графика напряжения на CONT с учётом внутреннего импеданса (25k min)
V_C = Vcc * (1 - exp(-t / tau_time_delay))
V_R = Vcc - V_C
I_R = V_R / R_time_delay
I_C = I_R
Z_C = simplify(V_C / I_C)
Z_Cont = 25e3;
Z_net = R_time_delay + 1 / (1/Z_C + 1/Z_Cont)
V_out = Vcc - R_time_delay * Vcc / Z_net


V_out_s = subs(V_out, {R_time_delay, C_time_delay, Vcc}, {R_time_delay_v, C_time_delay_v, Vcc_v})
figure("name", sprintf("Expected voltage at CONT pin, R = %g, C = %g", R_time_delay_v, C_time_delay_v));
ezplot(t, V_out_s)
xlim([0, R_time_delay_v * C_time_delay_v * time_delay_factor])
axis("square")
ylim([0, 12])

% Через примерно 3 секунды после включения напряжение на CONT станет около 10 вольт
t_fin = 3
V_Cont_fin = 10

% Компоненты, задающие режим работы таймера (astable mode)
Ra = 10e3; Rb = 1e3; C_pwm = 0.1e-6
freq_std = 1.44 / (C_pwm * (Ra + 2 * Rb))

% Частота и скважность при выбраном напряжении на CONT
tH = (Ra + Rb) * C_pwm * (log((Vcc - 1/2 * V_Cont_fin) / Vcc) - log((Vcc - V_Cont_fin) / Vcc))
tL = - Rb * C_pwm * log(1/2)

tH = double(subs(tH, {Vcc}, {Vcc_v}))

pwm_freq_at_v_cont = 1 / (tL + tH)
pwm_Q_at_v_cont = tH / (tL + tH)

% Желаемое напряжение и его изменение после интегрирующей RC-цепочки
V_target = 5; delta_V_target = 0.01; r_target = V_target / Vcc_v; k_target = (V_target + delta_V_target) / V_target;

k_expr = exp(tL / tau)
r_expr = (exp(tH / tau) - 1) / (exp((tH + tL) / tau) - 1)

tau_from_k = solve(k_target == k_expr, tau)
tau_from_k = double(tau_from_k)

% Фактическое напряжение и его дельта
r_actual = double(subs(r_expr, {tau}, {tau_from_k}))
k_actual = double(subs(k_expr, {tau}, {tau_from_k}))

V_actual = r_actual * Vcc_v

% Выбор параметров интегрирующей RC-цепочки и результриующее напряжение и его дельта после неё
R_int = 160E3; C_int = 0.22E-6;
tau_int_actual = R_int * C_int

r_actual = double(subs(r_expr, {tau}, {tau_int_actual}))
k_actual = double(subs(k_expr, {tau}, {tau_int_actual}))
V_actual = r_actual * Vcc_v
delta_V_actual = V_actual * k_actual - V_actual

% Выходной импеданс интегрирующей цепочки - для включения делителя напряжения после неё
Z_out_int = 1 / (1 / R_int + 1 / (1 / (2 * pi * pwm_freq_at_v_cont)))
Z_in_int = R_int + 1 / (2 * pi * pwm_freq_at_v_cont)

Time_graph = [0 : R_time_delay_v * C_time_delay_v * time_delay_factor / 100 : R_time_delay_v * C_time_delay_v * time_delay_factor];
Cont = double(subs(V_out_s, t, Time_graph));
%Cont = [0:0.1:Vcc_v];
tL_graph = -Rb * C_pwm * log(1/2);
tH_graph = (Ra + Rb) * C_pwm * (log((Vcc_v - 1/2 * Cont) / Vcc_v) - log((Vcc_v - Cont) / Vcc_v));
pwm_freq_graph = 1 ./ (tH_graph + tL_graph);
pwm_Q_graph = tH_graph ./ (tH_graph + tL_graph);
figure("name", sprintf("PWM freq, Ra = %g, Rb = %g, C = %g", Ra, Rb, C_pwm));
plot(Time_graph, pwm_freq_graph);
%set(gca, "xscale", "log");

figure("name", sprintf("PWM Q, Ra = %g, Rb = %g, C = %g", Ra, Rb, C_pwm));
plot(Time_graph, pwm_Q_graph);
%set(gca, "xscale", "log");
