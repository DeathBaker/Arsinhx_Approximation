#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <immintrin.h>
#include <math.h>

#define ln10 2.3025850929940456840179914546844
#define ln2 0.69314718055994530941723212145818

//  |number| <=2
double ln_small(double number)
{
    union
    {
        double num;
        int64_t bits;
    } temp;
    temp.num = number;
    int64_t exponent = ((temp.bits >> 52) & 0x7ff) - 0x3ff;
    temp.bits &= 0x800fffffffffffff;
    temp.bits |= 0x3ff0000000000000;
    double x = temp.num;
    // Polynomial of degree 7 from remez Approximant create with help of maple
    // source : https://de.maplesoft.com/support/help/maple/view.aspx?path=numapprox/remez
    // source :https://en.wikipedia.org/wiki/Minimax_approximation_algorithm
    temp.num = (-2.249637941 + (4.94490376300000000 + (-5.1945626890000005 + (4.007420810 + (-2.06908187600000005 + (0.6779732481 + (-0.1274994365 + 0.010484313120000005 * x) * x) * x) * x) * x) * x) * x) + exponent * ln2;
    return temp.num;
}
// newton Iteration Source: https://en.wikipedia.org/wiki/Newton%27s_method#Square_root
double newton_iteration(double number, double point, int count)
{
    if (count > 5)
    {
        return point;
    }
    point = point - (((point * point) - number) / (2 * point));
    return newton_iteration(number, point, count + 1);
}
// square root Source: https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Binary_estimates
double square_root(double number)
{
    union
    {
        double value;
        int64_t bits; /* data */
    } temp;
    temp.value = number;
    int64_t exponent = ((temp.bits >> 52) & 0x7ff) - 0x3ff;
    temp.bits &= 0x800fffffffffffff;
    if (exponent & 0x1)
    {
        temp.bits |= 0x4000000000000000;
    }
    else
    {
        temp.bits |= 0x3ff0000000000000;
    }
    exponent >>= 1;
    temp.value = newton_iteration(temp.value, 2, 0);
    int64_t final = ((((temp.bits >> 52) & 0x7ff) + exponent) << 52);
    temp.bits &= 0x800fffffffffffff;
    temp.bits |= final;
    return temp.value;
}
double approxArsinh_lookup(double x)
{
    if(x == INFINITY){
        return INFINITY;
    }

    else if (x == 0)
    {
        return 0.0;
    }
    else if(x == -INFINITY){
        return -INFINITY;
    }
    union
    {
        double value;
        int64_t bits;
    } temp;
    temp.value = x;
    int64_t exponent = ((temp.bits >> 52) & 0x7ff) - 0x3ff;
    if (exponent > 26)
    {
        return ln2 + ln_small(x);
    }
    else
    {
        return ln_small(x + square_root(x * x + 1));
    }
}

//  |number| <=2

double ln(double number)
{
    union
    {
        double data;
        int64_t bits;
    } temp;
    number--;
    temp.data = number;
    int64_t exponent = ((temp.bits >> 52) & 0x7ff) - 0x3ff;
    int count = exponent;
    double sum = 0;
    double sign = 1;
    double counter = 1;
    double x = number;
    while (count > -1022)
    {
        sum += (x / counter) * sign;
        x *= number;
        counter++;
        sign *= -1;
        count += exponent;
    }
    return sum;
}



double lessThanOne(double x, double exponent)
{
    double count = exponent;
    double power = x * x * x;
    double counter = 1;
    double nominator = 1;
    double denominator = 1;
    int sign = -1;
    double sum = x;
    while ((count < 1023) && counter < 86)
    {
        nominator *= (2 * counter - 1);
        denominator *= 2 * counter * (2 * counter + 1);
        sum += (sign) * ((nominator / denominator) * power);
        power *= x * x;
        counter++;
        sign *= -1;
        count += exponent;
    }
    return sum;
}
double arsinh_bigger_than_1(double number, double exponent)
{
    // our formel for ln is vaild only for x<=1 so we need to extract number of ln(10)
    double xLN = number;
    int power = 0;

    while (xLN >= 1)
    {
        xLN /= 10;
        power++;
    }

    // Setup all factors needed for ln calculation

    // setup other factors
    double sqr_number = -(number * number);
    int count = exponent;
    double sum = 0.0;
    double counter = 1.0;
    double dbl_fac_upper = 1;
    double dbl_fac_unter = 2;

    // loop if exponent of series is still >-53
    while (count < 1023 && counter < 151)
    {
        sum += (dbl_fac_upper / (2 * counter * dbl_fac_unter * (sqr_number)));
        counter++;
        dbl_fac_upper *= (2 * counter - 1);
        dbl_fac_unter *= (2 * counter);
        sqr_number *= -(number * number);
        count += exponent + 2;
    }

    //+ln(2) because we need ln(2|x|)
    return ln(xLN) + power * ln10 + ln2 - sum;
}

// ln(2|x|) (possibly not needed, see double arsinh_bigger_than_1(double number,int exponent) and double ln(double number) )



double approxArsinh_series(double x)
{
    if(x == INFINITY){
        return INFINITY;
    }
    else if(x == -INFINITY){
        return -INFINITY;
    }
    union
    {
        double p;
        int64_t bits;
    } temp;
    temp.p = x;
    int exponent = ((temp.bits >> 52) & 0x7ff) - 0x3ff;
    if (exponent < 0)
    {
        return lessThanOne(x, -exponent);
    }
    else if (x > 0)
    {
        return arsinh_bigger_than_1(x, exponent);
    }
    else
    {
        return -arsinh_bigger_than_1(-x, exponent);
    }
}

int main(){
    printf("%.18f",approxArsinh_series(50));
}



// using built-in arsinh Function



