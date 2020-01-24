# ref: http://www.math-cs.gordon.edu/courses/mat342/handouts/gauss.pdf
import numpy as np


class Matrix:
    def __init__(self, data, oldx, b):
        if len(data) == len(data[0]):
            self.data = data
            # Macierz musi być kwadratowa, czyli tu tylko parametr size.
            self.size = len(data)
            self.b = b
            self.oldx = oldx
            self.x = []
            for i in range(0, self.size):
                self.x.append(0)

    def __str__(self):
        temp = ""
        for y in self.data:
            for x in y:
                temp += " | " + str(x)
            temp += "\n"
        print("B = " + str(self.b))
        print("SOLVED X = " + str(self.x))
        print("GENERATED X = " + str(self.oldx))

        return temp

    def gauss_nopivoting(self):  # To dziala dobrze
        n = self.size
        data = self.data
        for k in range(0, n - 1):
            for i in range(k + 1, n):
                data[i, k] = data[i, k] / data[k, k]
                for j in range(k + 1, n):
                    data[i, j] = data[i, j] - (data[i, k] * data[k, j])
        for k in range(0, n - 1):
            for i in range(k + 1, n):
                self.b[i] = self.b[i] - (data[i, k] * self.b[k])
        for i in range(n - 1, -1, -1):
            s = self.b[i]
            for j in range(i + 1, n):
                s = s - (data[i, j] * self.x[j])
            self.x[i] = s / data[i, i]

    def gauss_partialpivoting(self):
        n = self.size
        a = self.data
        s = np.zeros(n)
        p = np.arange(n)
        for i in range(0, n):
            for j in range(0, n):
                temp = [s[i], abs(a[i, j])]
                s[i] = max(temp)
        for i in range(0, n - 1):
            rm = 0
            for j in range(i, n - 1):
                r = max(s[i], abs(a[p[j], i]))
                if r > rm:
                    rm = r
                    t = j
            temp = p[i]
            p[i] = p[t]
            p[t] = temp
            for j in range(i + 1, n):
                a[p[j], i] = -(a[p[j], i] / a[p[i], i])
                for k in range(i + 1, n):
                    a[p[j], k] += (a[p[j], i] * a[p[i], k])
        for i in range(0, n - 1):
            for j in range(i + 1, n):
                self.b[p[j]] += (a[p[j], i] * self.b[p[i]])
        for i in range(n - 1, -1, -1):
            s = self.b[p[i]]
            for j in range(i + 1, n):
                s -= (a[p[i], j] * self.x[j])
            self.x[i] = s / a[p[i], i]

    def gauss_completepivoting(self):
        n = self.size
        data = self.data
        col = np.arange(n)
        row = np.arange(n)
        for i in range(0, n):
            col[i] = i
            row[i] = i
        for k in range(0, n - 1):
            rm = 0
            for i in range(k, n):
                for j in range(k, n):
                    r = abs(data[row[i], col[j]])
                    if r > rm:
                        rm = r
                        rowindex = i
                        colindex = j
            temp = row[k]  # ZAMIANA INDEKSÓW WIERSZY.
            row[k] = row[rowindex]
            row[rowindex] = temp
            temp = col[k]  # ZAMIANA INDEKSÓW KOLUMN.
            col[k] = col[colindex]
            col[colindex] = temp
            for i in range(k + 1, n):
                data[row[i], col[k]] = data[row[i], col[k]] / data[row[k], col[k]]
                for j in range(k + 1, n):
                    data[row[i], col[j]] = data[row[i], col[j]] - (
                            data[row[i], col[k]] * data[row[k], col[j]])

        # forward elimination
        for k in range(0, n - 1):
            for i in range(k + 1, n):
                self.b[row[i]] = self.b[row[i]] - (data[row[i], col[k]] * self.b[row[k]])

        # solving
        for i in range(n - 1, -1, -1):
            s = self.b[row[i]]
            for j in range(i + 1, n):
                s = s - (data[row[i], col[j]] * self.x[col[j]])
            self.x[col[i]] = s / data[row[i], col[i]]

        # just in case
