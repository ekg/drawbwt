# drawbwt

## illustrating BWT based text search

Use cairo to build vector illustrations of BWT based searching on small strings. 

```
bin/drawbwt -o x.pdf -t mississippi -x 400 -y 400
```

To build, `cmake -H. -Bbuild && cmake --build build -- -j 4`.
