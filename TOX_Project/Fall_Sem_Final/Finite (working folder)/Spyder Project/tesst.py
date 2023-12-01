import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Generate some random data
data = np.random.randn(20)

# Set up the Seaborn line plot
sns.set()
sns.set_style("whitegrid")
ax = sns.lineplot(x=[0, 1, 2], y=[1, 2, 3])

# Set the ytick locations and labels, can also use np array here
ax.set_yticks([0, 1, 2, 3, 4])
ax.set_yticklabels(["A", "B", "C", "D", "E"])

for ind, label in enumerate(ax.get_xticklabels()):
    if ind % 5 == 0:  # every 10th label is kept
        label.set_visible(True)
    else:
        label.set_visible(False)

# Show the plot
plt.show()