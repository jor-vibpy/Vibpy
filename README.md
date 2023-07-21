# Vibpy

To install Vibpy, run the line below in a Terminal.

```
pip install -i https://test.pypi.org/simple/ Vibpy
```


## Utilizing the Package

This package allows the user to create Multi Degree of Freedom systems, obtaining the system characteristics and obtain the response of the system for a desiring input, allowing the plotting of the results.

### Examples

An example can be viwed in the code below.

```
  M = np.array([[1,0],[0,2]])
  K = np.array([[5,-2], [-2,3]])
  C = np.array([[4,-1], [-1,2]])
  
  CI = np.array([0.2,0,1,0])
  t = np.linspace(0,20,1001)
  F = np.array([1,2])
  w = np.array([3,3])
  
  sys = MDOF(M, K, C)
  sysn = sys.MDOF(CI, t, F, w, 'cos')
  
  sys.plot_responses(t)
```

The code above will return a plot with the displacement responses of the system.


For One Degree of Freedom systems, an example can be see below.

```
m = 0.0114
k = 430000
c = 1.400285685

CI = np.array([0.08,0])

t = np.linspace(0,1,1001)

F = 1

w = 2456.641553

sys2 = vp.DOF1(m, k, c)
sys3 = sys2.DOF1(CI, t, F, w)
sys2.plot_responses(t)
```
