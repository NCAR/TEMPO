import matplotlib.pyplot as plt


def get_data_from_file(fname):

    fn = open(fname, "r")
    first_line = fn.readline()
    var_names = list(first_line.split())

    list_of_lists = []
    for line in fn:
        split_line = list(line.strip().split())
        split_line_to_float = [float(e) for e in split_line]
        list_of_lists.append(split_line_to_float)

    transpose_list = list(zip(*list_of_lists))
    output_data = {}
    for i, var_name in enumerate(var_names):
        output_data[var_name] = transpose_list[i]
    return output_data


init_data = get_data_from_file("data/mpas_59lev_test.txt")
baseline_data = get_data_from_file("data/mpas_59lev_test_results_baseline.txt")
results_data = get_data_from_file("mpas_59lev_test_results.txt")

# plot
fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot([i*1000. for i in init_data['qr']], [i/100. for i in init_data['pressure']], c='k', label='init profile', linestyle='dashed', linewidth=0.5)
ax.plot([i*1000. for i in baseline_data['qr']], [i/100. for i in baseline_data['pressure']], c='k', label='baseline')
ax.plot([i*1000. for i in results_data['qr']], [i/100. for i in results_data['pressure']], c='grey', linestyle='dashed', label='update')
ax.set_ylim((1000, 100))
ax.legend()
ax.set_ylabel(r'Pressure [hPa]')
ax.set_xlabel(r'q$_{R}$ [g kg$^{-1}$]')
plt.savefig('mpas_59lev_test_results.pdf', bbox_inches='tight', dpi=25)

