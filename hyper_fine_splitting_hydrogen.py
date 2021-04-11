import numpy as np
import plotly.graph_objects as go

# General Parameters for both HFS and Lamb Shift
r_0 = 2.8179403262 * 10 ** -15  # meters
r_0 = 1  # atomic units
speed_of_light = 137.036  # atomic units

# HFS
delta_energy_0_hf = 1420.0  # MHz
delta_energy_hf = 1420.405751766710
# (10) MHz

delta_ratio = delta_energy_hf / delta_energy_0_hf
# 1.0002857406807817
a = 4 * (1 + delta_ratio)
# 8.001142962723126
b = -4

delta_0 = -(8 + 24 * delta_ratio)
# -32.00685777633876
delta_1 = 116 + (72 * delta_ratio) - 108 * (delta_ratio ** 2)
# 79.9588445240119
eulers_constant = 2 ** (-1 / 3) * (delta_1 + (delta_1 ** 2 - 4 * delta_0 ** 3) ** (1 / 2)) ** (1 / 3)
# 6.085965147441353
dt_electron_radius = -(r_0 / (3 * a)) * (b + eulers_constant + (delta_0 / eulers_constant))
# 3.7252086266951164e-16 # meters

# Round to the correct Precision of the NIST value for Hyper Fine Splitting
dt_electron_radius_sig_dig = float(np.format_float_scientific(dt_electron_radius, 10))
# 3.7252086267e-16 # meters
print(dt_electron_radius)

# Lamb Shift
g_p = 5.5856
c_prime = -(1 + r_0 / 2)
eulers_constant = 0.5772156649

b_prime = (1 + r_0 / 2 + r_0 ** 2 / 4)
a = (1 / 8 * np.pi) * (g_p / (2 * (3 ** 2) * speed_of_light ** 4))

predicted_lamb_shift = (a / 4) * (1 + 4 * c_prime + (r_0 ** 2 / 2 * dt_electron_radius_sig_dig ** 2) + (
        (4 * r_0 ** 2 - 3 * r_0) / dt_electron_radius_sig_dig)) + \
                       (a / 4) * ((r_0 ** 2 - 3 * r_0 - 4 * b_prime) * eulers_constant + (
        r_0 ** 2 - 3 * r_0 - 4 * b_prime) * np.log(dt_electron_radius_sig_dig))

print(predicted_lamb_shift)

predicted_lamb_shift_mhz = 4.13 * 10 ** 16 * predicted_lamb_shift * 1e-6

# Wikipedia says it should be 1047 MHz

measured_lamb_shift = 1047  # MHz

print(predicted_lamb_shift)

search_range = np.logspace(-4, 0, num=10000)
search_predicted_lamb_shift = []

for test_radius in search_range:
    print(test_radius)
    search_predicted_lamb_shift.append((a / 4) * (
            1 + 4 * c_prime + (r_0 ** 2 / 2 * test_radius ** 2) + ((4 * r_0 ** 2 - 3 * r_0) / test_radius)) + (
                                               a / 4) * ((r_0 ** 2 - 3 * r_0 - 4 * b_prime) * eulers_constant + (
            r_0 ** 2 - 3 * r_0 - 4 * b_prime) * np.log(test_radius)))

print(search_predicted_lamb_shift)
search_predicted_lamb_shift_mhz = np.array(search_predicted_lamb_shift) * 4.13 * 10 ** 16 * 1e-6

print(search_predicted_lamb_shift_mhz)
fig = go.Figure()

fig.add_trace(go.Scatter(x=search_range, y=search_predicted_lamb_shift_mhz, mode='markers', name='Search Radius'))

fig.add_hline(y=1047, line_dash="dot",
              annotation_text="Exerpimental Measurement",
              annotation_position="top right")

fig.add_trace(go.Scatter(x=[dt_electron_radius], y=[predicted_lamb_shift_mhz], name='Predicted Value', marker_symbol='star',
                         marker=dict(size=12,
                              line=dict(width=2,
                                        color='DarkSlateGrey'))))

fig.update_layout(title=
                  'Radius Search, Calculation, and Measurement Comparison',
                  xaxis_title='Radius of Electron (Atomic Units)', yaxis_title='Predicted Lamb Shift (MHz)')

fig.write_html('search.html')
