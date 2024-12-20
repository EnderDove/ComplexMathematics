﻿using System.Globalization;
using System.Numerics;

namespace ComplexMathematics
{
    /// <summary>
    /// A Class that provides Complex numbers and operations with them.
    /// Contains <paramref name="real"/> and <paramref name="imaginary"/> parts of number represented by double Type.
    /// </summary>
    [Serializable]
    public readonly struct Complex(double real, double imaginary) : IEquatable<Complex>, IFormattable, IAdditionOperators<Complex, Complex, Complex>,
        IAdditiveIdentity<Complex, Complex>, IDecrementOperators<Complex>, IDivisionOperators<Complex, Complex, Complex>, IEqualityOperators<Complex, Complex, bool>,
        IIncrementOperators<Complex>, IMultiplicativeIdentity<Complex, Complex>, IMultiplyOperators<Complex, Complex, Complex>, ISubtractionOperators<Complex, Complex, Complex>, IUnaryNegationOperators<Complex, Complex>
    {
        public double Real => _real;
        public double Imaginary => _imaginary;

        private readonly double _real = real;
        private readonly double _imaginary = imaginary;

        public double Argument => Math.Atan2(_imaginary, _real);
        public double Magnitude => Abs(this);
        public Complex Conjugated => Conjugate(this);

        #region Constants
        public static Complex Zero => new(0, 0);
        public static Complex AdditiveIdentity => Zero;
        public static Complex RealOne => new(1, 0);
        public static Complex MultiplicativeIdentity => RealOne;
        public static Complex ImaginaryOne => new(0, 1);
        public static Complex EpsilonForAnalysis => new(0.00001d, 0);
        public static Complex EulersConstant => new(0.5772157, 0);
        public static Complex PI => new(Math.PI, 0);
        public static Complex Tau => new(Math.Tau, 0);
        public static Complex E => new(Math.E, 0);
        #endregion

        #region Generators
        public static Complex FromPolarCoordinates(double magnitude, double argument) => new(magnitude * Math.Cos(argument), magnitude * Math.Sin(argument));

        public static implicit operator Complex(short value) => new(value, 0);

        public static implicit operator Complex(int value) => new(value, 0);

        public static implicit operator Complex(long value) => new(value, 0);

        public static implicit operator Complex(ushort value) => new(value, 0);

        public static implicit operator Complex(uint value) => new(value, 0);

        public static implicit operator Complex(ulong value) => new(value, 0);

        public static implicit operator Complex(sbyte value) => new(value, 0);

        public static implicit operator Complex(byte value) => new(value, 0);

        public static implicit operator Complex(float value) => new(value, 0);

        public static implicit operator Complex(double value) => new(value, 0);

        public static explicit operator Complex(BigInteger value) => new((double)value, 0);

        public static explicit operator Complex(decimal value) => new((double)value, 0);

        #endregion

        #region ComplexFunctions
        public static double Abs(Complex value)
        {
            if (double.IsInfinity(value._real) || double.IsInfinity(value._imaginary))
                return double.PositiveInfinity;

            var real = Math.Abs(value._real);
            var imaginary = Math.Abs(value._imaginary);
            if (real > imaginary)
            {
                var ratio1 = imaginary / real;
                return real * Math.Sqrt(1 + ratio1 * ratio1);
            }

            if (imaginary == 0)
                return real;

            var ratio2 = real / imaginary;
            return imaginary * Math.Sqrt(1 + ratio2 * ratio2);
        }

        public static Complex Conjugate(Complex value)
        {
            return new Complex(value._real, -value._imaginary);
        }

        public static Complex Log(Complex value)
        {
            return new(Math.Log(value.Magnitude), value.Argument);
        }

        public static Complex Log(Complex value, Complex baseValue)
        {
            return Log(value) / Log(baseValue);
        }

        public static Complex Exp(Complex value)
        {
            return FromPolarCoordinates(Math.Exp(value._real), value._imaginary);
        }

        public static Complex Pow(Complex value, Complex power)
        {
            if (power == Zero)
                return RealOne;

            if (value == Zero)
                return Zero;

            var valueArgument = value.Argument;
            var valueMagnitude = value.Magnitude;
            var resultArgument = power._real * valueArgument + power._imaginary * Math.Log(valueMagnitude);
            var resultMagnitude = Math.Pow(valueMagnitude, power._real) * Math.Exp(-power._imaginary * valueArgument);
            return FromPolarCoordinates(resultMagnitude, resultArgument);
        }

        public static Complex Sqrt(Complex value)
        {
            return FromPolarCoordinates(Math.Sqrt(value.Magnitude), value.Argument / 2);
        }

        public static Complex Sin(Complex value)
        {
            return new(Math.Sin(value._real) * Math.Cosh(value._imaginary), Math.Cos(value._real) * Math.Sinh(value._imaginary));
        }

        public static Complex Sinh(Complex value)
        {
            return new(Math.Sinh(value._real) * Math.Cos(value._imaginary), Math.Cosh(value._real) * Math.Sin(value._imaginary));
        }

        public static Complex Asin(Complex value)
        {
            return -ImaginaryOne * Log(ImaginaryOne * value + Sqrt(RealOne - value * value));
        }

        public static Complex Asinh(Complex value)
        {
            return Log(value + Sqrt(RealOne + value * value));
        }

        public static Complex Cos(Complex value)
        {
            return new(Math.Cos(value._real) * Math.Cosh(value._imaginary), -Math.Sin(value._real) * Math.Sinh(value._imaginary));
        }

        public static Complex Cosh(Complex value)
        {
            return new(Math.Cosh(value._real) * Math.Cos(value._imaginary), Math.Sinh(value._real) * Math.Sin(value._imaginary));
        }

        public static Complex Acos(Complex value)
        {
            return -ImaginaryOne * Log(value + ImaginaryOne * Sqrt(RealOne - value * value));
        }

        public static Complex Acosh(Complex value)
        {
            return Log(ImaginaryOne * (-value + Sqrt(RealOne + value * value)));
        }

        public static Complex Tan(Complex value)
        {
            return Sin(value) / Cos(value);
        }

        public static Complex Tanh(Complex value)
        {
            return Sinh(value) / Cosh(value);
        }

        public static Complex Atan(Complex value)
        {
            return 0.5 * ImaginaryOne * (Log(RealOne - ImaginaryOne * value) - Log(RealOne + ImaginaryOne * value));
        }

        public static Complex Atanh(Complex value)
        {
            return -0.5 * (Log(RealOne - value) - Log(RealOne + value));
        }

        public static Complex Lerp(Complex first, Complex second, double factor)
        {
            return first + (second - first) * factor;
        }

        public static double InverseLerp(Complex first, Complex second, Complex value)
        {
            return ((first - value) / (first - second)).Real;
        }

        public static Complex SLerp(Complex first, Complex second, double factor)
        {
            var ratio = second / first;
            return first * FromPolarCoordinates(1 + (ratio.Magnitude - 1) * factor, ratio.Argument * factor);
        }

        public static Complex Gamma(Complex value)
        {
            var result = value * Exp(EulersConstant * value);
            for (int k = 1; k < 1000; k++)
                result *= (1 + value / k) * Exp(-value / k);
            return 1 / result;
        }

        public static Complex Beta(Complex a, Complex b)
        {
            return Gamma(a) * Gamma(b) / Gamma(a + b);
        }

        public static Complex Zeta(Complex value)
        {
            if (value == Zero) return -0.5;
            if (value == RealOne) return new(double.PositiveInfinity, 0);
            if (value._real < 0.5)
                return Pow(2, value) * Pow(PI, value - 1) * Sin(0.5 * Math.PI * value) * Gamma(1 - value) * Zeta(1 - value);
            var sum = Zero;
            for (int k = 1; k < 1000; k++)
                sum += k * (k + 1) / 2 * ((2 * k + 3 + value) / Pow(k + 1, value + 2) - (2 * k - 1 - value) / Pow(k, value + 2));
            return sum / (value - 1);
        }

        public static Complex LambertW(Complex value)
        {
            var w = Zero;
            for (int k = 0; k < 1000; k++)
            {
                var wTimesExpW = w * Exp(w);
                w -= (wTimesExpW - value) / ((w + 1) * Exp(w) - (w + 2) * (wTimesExpW - value) / (2 * w + 2));
            }
            return w;
        }
        #endregion

        #region Transforms
        /// <summary>
        /// Fast Fourier Transform Algorithm with 0-padding for not 2-th power inputs.
        /// </summary>
        /// <returns> Complex[] with length of input array scaled to minimal power of 2 greater then k </returns>
        public static Complex[] FastFourierTransform(Complex[] input)
        {
            var initialPointsCount = input.Length;
            if (initialPointsCount <= 1) return input;
            var pointsCount = (int)Math.Pow(2, (int)Math.Log2(initialPointsCount - 1) + 1);

            var even = new Complex[pointsCount / 2];
            var odd = new Complex[pointsCount / 2];
            for (int i = 0; i < initialPointsCount / 2; i++)
            {
                even[i] = input[2 * i];
                odd[i] = input[2 * i + 1];
            }

            even = FastFourierTransform(even);
            odd = FastFourierTransform(odd);
            var result = new Complex[pointsCount];
            for (int frequency = 0; frequency < pointsCount / 2; frequency++)
            {
                var value = FromPolarCoordinates(1, -Math.Tau * frequency / pointsCount) * odd[frequency];
                result[frequency] = even[frequency] + value;
                result[frequency + pointsCount / 2] = even[frequency] - value;
            }
            return result;
        }

        /// <summary>
        /// Discrete Fourier Transform Algorithm. Always prefer to use FastFourierTransform.
        /// </summary>
        /// <returns> Complex[] with length of input array</returns>
        public static Complex[] DiscreteFourierTransform(Complex[] input)
        {
            int pointsCount = input.Length;
            var result = new Complex[pointsCount];
            for (int point = 0; point < pointsCount; point++)
            {
                Complex coefficient = 0;
                for (int frequency = 0; frequency < pointsCount; frequency++)
                    coefficient += input[frequency] * FromPolarCoordinates(1, -Math.Tau * frequency * point / pointsCount);
                result[point] = coefficient;
            }
            return result;
        }

        /// <summary>
        /// Inverse Fast Fourier Transform Algorithm with 0-padding for not 2-th power inputs.
        /// </summary>
        /// <returns> Complex[] with length of input array scaled to minimal power of 2 greater then k </returns>
        public static Complex[] InverseFastFourierTransform(Complex[] input)
        {
            var pointsCount = input.Length;
            var result = new Complex[pointsCount];
            for (int i = 0; i < pointsCount; i++)
                result[i] = result[i].Conjugated;
            result = FastFourierTransform(result);
            for (int i = 0; i < pointsCount; i++)
                result[i] = result[i].Conjugated;
            return result;
        }

        /// <summary>
        /// Inverse Discrete Fourier Transform Algorithm. Always prefer to use InverseFastFourierTransform.
        /// </summary>
        /// <returns> Complex[] with length of input array</returns>
        public static Complex[] InverseDiscreteFourierTransform(Complex[] input)
        {
            var pointsCount = input.Length;
            var result = new Complex[pointsCount];
            for (int i = 0; i < pointsCount; i++)
                result[i] = result[i].Conjugated;
            result = DiscreteFourierTransform(result);
            for (int i = 0; i < pointsCount; i++)
                result[i] = result[i].Conjugated;
            return result;
        }
        #endregion

        #region ToString
        public override string ToString()
        {
            return string.Format(CultureInfo.CurrentCulture, $"{_real} + {_imaginary}i");
        }

        public string ToString(string format)
        {
            var real = this._real.ToString(format, CultureInfo.CurrentCulture);
            var imaginary = this._imaginary.ToString(format, CultureInfo.CurrentCulture);
            return string.Format(CultureInfo.CurrentCulture, $"{real} + {imaginary}i");
        }

        public string ToString(IFormatProvider provider)
        {
            return string.Format(provider, $"{_real} + {_imaginary}i");
        }

        public string ToString(string? format, IFormatProvider? provider)
        {
            var real = this._real.ToString(format, CultureInfo.CurrentCulture);
            var imaginary = this._imaginary.ToString(format, CultureInfo.CurrentCulture);
            return string.Format(provider, $"{real} + {imaginary}i");
        }
        #endregion

        #region Operators
        public static Complex operator +(Complex left, Complex right)
        {
            return new Complex(left._real + right._real, left._imaginary + right._imaginary);
        }

        public static Complex operator ++(Complex value)
        {
            return new(value._real + 1, value._imaginary);
        }

        public static Complex operator -(Complex value)
        {
            return new(-value._real, -value._imaginary);
        }

        public static Complex operator -(Complex left, Complex right)
        {
            return new(left._real - right._real, left._imaginary - right._imaginary);
        }

        public static Complex operator --(Complex value)
        {
            return new(value._real - 1, value._imaginary);
        }

        public static Complex operator *(Complex left, Complex right)
        {
            var real = left._real * right._real - left._imaginary * right._imaginary;
            var imaginary = left._imaginary * right._real + left._real * right._imaginary;
            return new(real, imaginary);
        }

        public static Complex operator *(Complex left, double right)
        {
            return new(left._real * right, left._imaginary * right);
        }

        public static Complex operator *(double left, Complex right)
        {
            return new(left * right._real, left * right._real);
        }

        public static Complex operator /(Complex left, Complex right)
        {
            if (Math.Abs(right._imaginary) < Math.Abs(right._real))
            {
                var num = right._imaginary / right._real;
                return new((left._real + left._imaginary * num) / (right._real + right._imaginary * num), (left._imaginary - left._real * num) / (right._real + right._imaginary * num));
            }
            var num2 = right._real / right._imaginary;
            return new((left._imaginary + left._real * num2) / (right._imaginary + right._real * num2), (-left._real + left._imaginary * num2) / (right._imaginary + right._real * num2));
        }

        public static Complex operator /(Complex left, double right)
        {
            return new(left._real / right, left._imaginary / right);
        }

        public static bool operator ==(Complex left, Complex right)
        {
            return left._real == right._real && left._imaginary == right._imaginary;
        }

        public static bool operator ==(double left, Complex right)
        {
            return left == right._real && right._imaginary == 0;
        }

        public static bool operator ==(Complex left, double right)
        {
            return left._real == right && left._imaginary == 0;
        }

        public static bool operator !=(Complex left, Complex right)
        {
            return left._real != right._real || left._imaginary != right._imaginary;
        }

        public static bool operator !=(double left, Complex right)
        {
            return left != right._real || right._imaginary != 0;
        }

        public static bool operator !=(Complex left, double right)
        {
            return left._real != right || left._imaginary != 0;
        }
        #endregion

        #region Equals
        public override bool Equals(object? obj)
        {
            if (obj is not Complex)
                return false;
            return this == (Complex)obj;
        }

        public bool Equals(Complex value)
        {
            return _real.Equals(value._real) && _imaginary.Equals(value._imaginary);
        }
        public override int GetHashCode()
        {
            int num = 99999997;
            int num2 = _real.GetHashCode() % num;
            int hashCode = _imaginary.GetHashCode();
            return num2 ^ hashCode;
        }
        #endregion
    }
}
