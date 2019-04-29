# drawbwt

## illustrating BWT based text search

Use cairo to build vector illustrations of BWT based searching on small strings. 

```
bin/drawbwt -o x.pdf -t mississippi -x 350 -y 930 -s iss
```

![searching for iss in mississippi](https://raw.githubusercontent.com/ekg/drawbwt/master/img/mississippi.png)

To build, `cmake -H. -Bbuild && cmake --build build -- -j 4`.
