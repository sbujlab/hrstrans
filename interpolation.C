double interpolation(double x, double x0, double y0, double x1, double y1){

double slope = y1-y0/x1-x0;


double y = slope*(x-x0)+y0;



return y;

}

double loginterpolation(double x, double x0, double y0, double x1, double y1){

double slope = (log(x) - log(x0))/(log(x1)-log(x0));

double y = slope*(x1-x0)+y0;

return y;


}

