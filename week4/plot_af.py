import matplotlib.pyplot as plt
import pandas as pd
import sys

def hist_af(df):
    data = df.iloc[:,4:]
    names = data.columns

    data.hist()
    plt.suptitle('AF dist')
    plt.show()

if __name__ == "__main__":
    file = sys.argv[1]
    df = pd.read_csv(file,sep='\t')
    hist_af(df)
