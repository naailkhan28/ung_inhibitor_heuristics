import matplotlib.pyplot as plt

hit_rates = [0, 4.94, 6.14, 2.21, 2.68, 0.97]

labels = ["<25%", "25 - 30%", "30 - 35%", "35 - 40%", "40 - 45%", "45%<"]

chart = plt.bar(labels, hit_rates, color="#91BBC2", edgecolor="black")

plt.ylabel("Hit Rate [%]")
plt.xlabel("GC Content [%]")

plt.show()