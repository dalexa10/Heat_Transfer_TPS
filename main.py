import pickle
import matplotlib.ticker as mtick
from heat_transfer import transient_heat_transfer
from constructors.materials import Multi_Layer_TPS
import matplotlib.pyplot as plt
from constructors.data.TPS_layouts_data import TPS_dict_B
from constructors.fluid_states import Fluid_State_CPG

TPS_instance = Multi_Layer_TPS(TPS_dict_B)

# ------------------------------------------------------------
# Load flow output (object) obtained from combustor
# ------------------------------------------------------------
with open('flow_data/flow_data.pkl', 'rb') as inp:
    f_cb_dict = pickle.load(inp)

f_cb = Fluid_State_CPG(**f_cb_dict)

T_v_ls, q_conv_ls, t_ls, tot_time = transient_heat_transfer(TPS=TPS_instance,
                                                            flow=f_cb,
                                                            Tw0=300,
                                                            tf=3,
                                                            dt=0.05,
                                                            reg='lam')
print("Wall temperature is {:.2f} at time {:.1f}".format(T_v_ls[-1][0][0], tot_time))

fig, ax = plt.subplots(1, 2)
ax[0].plot(T_v_ls[10], label='t = {:.2f} [s]'.format(t_ls[10]))
ax[0].plot(T_v_ls[20], label='t = {:.2f} [s]'.format(t_ls[20]))
ax[0].plot(T_v_ls[30], label='t = {:.2f} [s]'.format(t_ls[30]))
ax[0].plot(T_v_ls[-1], label='t = {:.2f} [s]'.format(t_ls[-1]))
ax[0].set_ylabel('Temperature [K]')
ax[0].set_xlabel('Nodes')
ax[0].legend(loc='best')

ax[1].plot(t_ls, q_conv_ls, '-k')
ax[1].set_xlabel('time [s]')
ax[1].set_ylabel('Heat transfer [W/$m^2$]')
formatter = mtick.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 1))
ax[1].yaxis.set_major_formatter(formatter)
plt.tight_layout()
plt.show()

