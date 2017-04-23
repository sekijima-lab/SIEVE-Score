import matplotlib
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pandas as pd


df = pd.read_csv("0309_fiximportance.importance",
                    names=["name", "value"], index_col=0)
df = df.sort_values(by=["value"], ascending=False)[:10]

plt.figure()
ax = df.plot(kind="bar", alpha=0.75, rot=30, legend=False, fontsize=14)
ax.set_xlabel("")
plt.savefig("importance.png",bbox_inches='tight')
