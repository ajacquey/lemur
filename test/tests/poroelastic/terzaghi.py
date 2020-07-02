import os, csv
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

def analyticalSolution(z, t):
  res = 0.0
  for k in range(1,500):
    res += 4.0 / np.pi * np.power(-1.0, k - 1.0) / (2.0 * k - 1.0) * np.cos((2.0 * k -1.0) * np.pi / 2.0 * z) * np.exp(-np.power(2.0 * k - 1.0, 2.0) * np.power(np.pi, 2.0) / 4.0 * t)
  return res

def parseCsvFile(filename, column_keys=None):
    column_index = {} # mapping, key=column_key, value=corresponding column index
    data = {} # dict of data, key=column_key, value=data list (floats)
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        line_i = 0 # line index
        for row in csvreader:
            if line_i == 0:
                # Headers. Find interesting columns
                headers = row
                # prepare structure for all columns we want
                if column_keys is None:
                    # grab all data in file
                    column_keys_we_want = [elt.lower() for elt in headers]
                else:
                    # grab only requested data from file
                    assert type(column_keys)==type([])
                    column_keys_we_want = column_keys
                for column_key in column_keys_we_want:
                    data[column_key] = []
                for column_i, elt in enumerate(headers):
                    elt_lower = elt.lower()
                    if elt_lower in column_keys_we_want:
                        column_index[elt_lower] = column_i
                line_i += 1
                continue
            # Data line
            if len(row) < len(headers):
                break # finished reading all data
            for column_key in column_keys_we_want:
                data[column_key].append(float(row[column_index[column_key]]))
            line_i += 1
            continue # go to next data line
    return data

def numericalSolution(t):
  if t==0.001:
    filename = 'terzaghi_csv_line_pf_0010.csv'
  elif t==0.01:
    filename = 'terzaghi_csv_line_pf_0100.csv'
  elif t==0.05:
    filename = 'terzaghi_csv_line_pf_0141.csv'
  elif t==0.1:
    filename = 'terzaghi_csv_line_pf_0191.csv'
  elif t==0.5:
    filename = 'terzaghi_csv_line_pf_0232.csv'
  elif t==1.0:
    filename = 'terzaghi_csv_line_pf_0282.csv'
  else:
    print('Unknown filename!')
    exit

  data = parseCsvFile(filename)
  z = np.asarray(data['z'])
  p = np.divide(np.asarray(data['pf']), p0)

  return (z, p)
    
if __name__ == "__main__":
  h = 10.0
  phi = 0.1
  k = 1.5
  mu = 1.0
  Kf = 8.0
  Ks = 10.0
  K = 4.0
  G = 3.0
  alpha = 1.0 - K / Ks
  S = phi / Kf + (alpha - phi) / Ks
  mv = 1.0 / (K + 4.0 / 3.0 * G)
  Cv = k / (mu * (S + alpha**2 * mv))
  p0 = alpha * mv / (S + alpha**2 *mv)

  z = np.linspace(0.0, 1.0, 100)

  # Line 1 Cv * t / h**2 = 0.001
  l1 = analyticalSolution(z, 0.001)
  (zn1, ln1) = numericalSolution(0.001)

  # Line 2 Cv * t / h**2 = 0.01
  l2 = analyticalSolution(z, 0.01)
  (zn2, ln2) = numericalSolution(0.01)

  # Line 3 Cv * t / h**2 = 0.05
  l3 = analyticalSolution(z, 0.05)
  (zn3, ln3) = numericalSolution(0.05)

  # Line 4 Cv * t / h**2 = 0.1
  l4 = analyticalSolution(z, 0.1)
  (zn4, ln4) = numericalSolution(0.1)

  # Line 5 Cv * t / h**2 = 0.5
  l5 = analyticalSolution(z, 0.5)
  (zn5, ln5) = numericalSolution(0.5)

  # Line 6 Cv * t / h**2 = 1.0
  l6 = analyticalSolution(z, 1.0)
  (zn6, ln6) = numericalSolution(1.0)

  fig, ax = plt.subplots()
  fig.set_size_inches(6, 5)
  ax.plot(l1, z, lw=1, color='xkcd:blue', label='0.001')
  ax.plot(l2, z, lw=1, color='xkcd:blue', label='0.01')
  ax.plot(l3, z, lw=1, color='xkcd:blue', label='0.05')
  ax.plot(l4, z, lw=1, color='xkcd:blue', label='0.1')
  ax.plot(l5, z, lw=1, color='xkcd:blue', label='0.5')
  ax.plot(l6, z, lw=1, color='xkcd:blue', label='1.0')
  ax.plot(ln1, zn1, ls='none', marker='o', ms=4, color='k')
  ax.plot(ln2, zn2, ls='none', marker='o', ms=4, color='k')
  ax.plot(ln3, zn3, ls='none', marker='o', ms=4, color='k')
  ax.plot(ln4, zn4, ls='none', marker='o', ms=4, color='k')
  ax.plot(ln5, zn5, ls='none', marker='o', ms=4, color='k')
  ax.plot(ln6, zn6, ls='none', marker='o', ms=4, color='k')
  # ax.set_xlim(0.0, 1.0)
  # ax.set_ylim(0.0, 1.0)
  ax.grid()
  ax.set_xlabel(r'$\frac{p_{f}}{p_{0}}$', fontsize=16)
  ax.set_ylabel(r'$\frac{z}{h}$', rotation=0, fontsize=16)
  ax.set_title('Terzaghi\'s consolidation problem')

  plt.savefig('terzaghi.png', dpi=200, bbox_inches='tight')