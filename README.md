# drawbwt

## illustrating BWT based text search

Vector illustrations of BWT based searching on small strings, configurable for arbitrary text and search patterns.
The search process rendering is inspired by [Alex Bowe's representation of searching in the FM-index](https://alexbowe.com/fm-index/).

Usage:

```
git clone https://github.com/ekg/drawbwt.git
cd drawbwt
cmake -H. -Bbuild && cmake --build build -- -j 4
bin/drawbwt -o x.pdf -t mississippi -x 350 -y 930 -s iss
```

This yields the following vector graphic (here rendered in PNG for the web):

![searching for iss in mississippi](https://raw.githubusercontent.com/ekg/drawbwt/master/img/mississippi.png)

