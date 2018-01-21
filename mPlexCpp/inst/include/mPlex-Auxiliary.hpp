

/*
expand.grid should end up looking something like this in c++

template<class T iter>
void expand_grid(T iterOut, T iterBegin, T iterEnd){
  do the stuff ...
};

template<class T, class... Args>
void expand_grid(T iterOut, T iterBegin, T iterEnd, Args... args) {
  expand_grid(iterOut, iterBegin, iterEnd);
  expand_grid(args);
}
*/
