import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def plot_ptm(ptm, title = ''):
    matrix = np.real_if_close(ptm.full())
    if matrix.shape[0] == 4:
        label = ['I','X','Y','Z']
        fig, (f1) = plt.subplots(
            nrows = 1, ncols = 1,
            figsize=(5, 4))
    elif matrix.shape[0] == 16:
        label = ['II','IX','IY','IZ','XI','XX','XY','XZ','YI','YX','YY','YZ','ZI','ZX','ZY','ZZ']
        fig, (f1) = plt.subplots(
            nrows = 1, ncols = 1,
            figsize=(6, 5))
    else:
        print('ERROR in size of ptm')

    fontsizes = 12
    cmap_set = 'RdBu'#'bwr'#'RdBu'
    cb = f1.imshow(matrix, cmap = cmap_set,vmax = 1,vmin = -1)
    fig.colorbar(cb, ax=f1, ticks=[-1, -1/2, 0, 1/2, 1])

    f1.xaxis.set_major_locator(ticker.MultipleLocator(1))
    f1.yaxis.set_major_locator(ticker.MultipleLocator(1))


    f1.set_xticklabels([''] + label, rotation=45, fontsize=fontsizes)
    f1.set_yticklabels([''] + label, fontsize=fontsizes)
    f1.set_title(title, fontsize=np.around(1.3*fontsizes))
    plt.show()
    
    
    
def plot_ptm_compare(ptm1, ptm2, dif = 'False', title1 = '', title2 = '', title3 = ''):
    matrix1 = np.real_if_close(ptm1.full())
    matrix2 = np.real_if_close(ptm2.full())
    figsize_coef = 1
    if dif == "True":
        figsize_coef = 1.5
        if matrix1.shape[0] == 4:
            label = ['I','X','Y','Z']
            fig, (f1, f2, f3) = plt.subplots(
                nrows = 1, ncols = 3,
                figsize=(12*figsize_coef, 4))
        elif matrix1.shape[0] == 16:
            label = ['II','IX','IY','IZ','XI','XX','XY','XZ','YI','YX','YY','YZ','ZI','ZX','ZY','ZZ']
            fig, (f1, f2, f3) = plt.subplots(
                nrows = 1, ncols = 3,
                figsize=(14*figsize_coef, 5))
        else:
            print('ERROR in size of ptm')
    else:
        if matrix1.shape[0] == 4:
            label = ['I','X','Y','Z']
            fig, (f1, f2) = plt.subplots(
                nrows = 1, ncols = 2,
                figsize=(12*figsize_coef, 4))
        elif matrix1.shape[0] == 16:
            label = ['II','IX','IY','IZ','XI','XX','XY','XZ','YI','YX','YY','YZ','ZI','ZX','ZY','ZZ']
            fig, (f1, f2) = plt.subplots(
                nrows = 1, ncols = 2,
                figsize=(14*figsize_coef, 5))
        elif matrix1.shape[0] == 64:
            label = ['III', 'IIX', 'IIY', 'IIZ', 'IXI', 'IXX', 'IXY', 'IXZ', 'IYI', 'IYX', 'IYY', 'IYZ', 'IZI', 'IZX', 'IZY', 'IZZ',
                     'XII', 'XIX', 'XIY', 'XIZ', 'XXI', 'XXX', 'XXY', 'XXZ', 'XYI', 'XYX', 'XYY', 'XYZ', 'XZI', 'XZX', 'XZY', 'XZZ',
                     'YII', 'YIX', 'YIY', 'YIZ', 'YXI', 'YXX', 'YXY', 'YXZ', 'YYI', 'YYX', 'YYY', 'YYZ', 'YZI', 'YZX', 'YZY', 'YZZ',
                     'ZII', 'ZIX', 'ZIY', 'ZIZ', 'ZXI', 'ZXX', 'ZXY', 'ZXZ', 'ZYI', 'ZYX', 'ZYY', 'ZYZ', 'ZZI', 'ZZX', 'ZZY', 'ZZZ',
                     ]
            fig, (f1, f2) = plt.subplots(
                nrows=1, ncols=2,
                figsize=(14 * figsize_coef, 5))
        else:
            print('ERROR in size of ptm')

    fontsizes = 12
    cmap_set = 'RdBu'#'bwr'#'RdBu'
    cb = f1.imshow(matrix1, cmap = cmap_set,vmax = 1,vmin = -1)
    fig.colorbar(cb, ax=f1, ticks=[-1, -1/2, 0, 1/2, 1])
    f1.xaxis.set_major_locator(ticker.MultipleLocator(1))
    f1.yaxis.set_major_locator(ticker.MultipleLocator(1))
    f1.set_xticklabels([''] + label, rotation=45, fontsize=fontsizes)
    f1.set_yticklabels([''] + label, fontsize=fontsizes)
    f1.set_title(title1, fontsize=np.around(1.3*fontsizes))
    
    cb = f2.imshow(matrix2, cmap = cmap_set,vmax = 1,vmin = -1)
    fig.colorbar(cb, ax=f2, ticks=[-1, -1/2, 0, 1/2, 1])
    f2.xaxis.set_major_locator(ticker.MultipleLocator(1))
    f2.yaxis.set_major_locator(ticker.MultipleLocator(1))
    f2.set_xticklabels([''] + label, rotation=45, fontsize=fontsizes)
    f2.set_yticklabels([''] + label, fontsize=fontsizes)
    f2.set_title(title2, fontsize=np.around(1.3*fontsizes))
    
    if dif == "True":
        m = matrix2-matrix1
        cb = f3.imshow(m, cmap = cmap_set,vmax = np.max(np.abs(m)),vmin = -np.max(np.abs(m)))
        fig.colorbar(cb, ax=f3)
        f3.xaxis.set_major_locator(ticker.MultipleLocator(1))
        f3.yaxis.set_major_locator(ticker.MultipleLocator(1))
        f3.set_xticklabels([''] + label, rotation=45, fontsize=fontsizes)
        f3.set_yticklabels([''] + label, fontsize=fontsizes)
        f3.set_title(title3, fontsize=np.around(1.3*fontsizes))
    plt.show()

    d = np.sqrt(matrix1.shape[0])
    print('Fidelity', (np.trace(matrix1@matrix2.T) + d)/(d*(d+1)))