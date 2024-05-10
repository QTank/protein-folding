from EnergyTable import read_data

# energy_table, index_dict = read_data()

# print(energy_table[index_dict['Y']][index_dict['T']])
# print(energy_table[index_dict['P']][index_dict['Y']])
# print(energy_table[index_dict['P']][index_dict['G']])

a = []
for i in range(8):
    for j in range(8 - i):
        print(f"{i}, {j}, {i**2 + j**2}")
        a.append(i**2 + j**2)

print(a.sort())
a1 = []

for i in a:
    if i not in a1:
        a1.append(i)

print(a1)


