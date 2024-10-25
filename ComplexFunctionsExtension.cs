namespace ComplexMathematics
{
    /// <summary>
    /// A class that provides extension methods for complex-valued functions.
    /// </summary>
    public static class ComplexFunctionsExtension
    {
        #region Partial Derivatives
        /// <summary>
        /// Partial Derivative of U with respect to x where f(z) = U(x, y) + i*V(x, y)
        /// </summary>
        public static double PartialDerivativeUx(this Func<Complex, Complex> function, Complex point)
        {
            var endPoint = new Complex(point.Real + Complex.EpsilonForAnalysis.Real, point.Imaginary);
            return (function(endPoint).Real - function(point).Real) / Complex.EpsilonForAnalysis.Real;
        }

        /// <summary>
        /// Partial Derivative of U with respect to y where f(z) = U(x, y) + i*V(x, y)
        /// </summary>
        public static double PartialDerivativeUy(this Func<Complex, Complex> function, Complex point)
        {
            var endPoint = new Complex(point.Real, point.Imaginary + Complex.EpsilonForAnalysis.Real);
            return (function(endPoint).Real - function(point).Real) / Complex.EpsilonForAnalysis.Real;
        }

        /// <summary>
        /// Partial Derivative of V with respect to x where f(z) = U(x, y) + i*V(x, y)
        /// </summary>
        public static double PartialDerivativeVx(this Func<Complex, Complex> function, Complex point)
        {
            var endPoint = new Complex(point.Real + Complex.EpsilonForAnalysis.Real, point.Imaginary);
            return (function(endPoint).Imaginary - function(point).Imaginary) / Complex.EpsilonForAnalysis.Real;
        }

        /// <summary>
        /// Partial Derivative of V with respect to y where f(z) = U(x, y) + i*V(x, y)
        /// </summary>
        public static double PartialDerivativeVy(this Func<Complex, Complex> function, Complex point)
        {
            var endPoint = new Complex(point.Real, point.Imaginary + Complex.EpsilonForAnalysis.Real);
            return (function(endPoint).Imaginary - function(point).Imaginary) / Complex.EpsilonForAnalysis.Real;
        }
        #endregion

        #region Whole Derivatives
        public static bool IsDifferentiable(this Func<Complex, Complex> function, Complex point)
        {
            return (Math.Abs(function.PartialDerivativeUx(point) - function.PartialDerivativeVy(point)) <= Complex.EpsilonForAnalysis.Real * 2) &&
                (Math.Abs(function.PartialDerivativeUy(point) + function.PartialDerivativeVx(point)) <= Complex.EpsilonForAnalysis.Real * 2);
        }

        /// <summary>
        /// Complex Derivative df(z) = (U'x + i*V'x)dz where f(z) = U(x, y) + i*V(x, y)
        /// </summary>
        public static Complex Derivative(this Func<Complex, Complex> function, Complex point)
        {
            return new(function.PartialDerivativeUx(point), function.PartialDerivativeVx(point));
        }
        #endregion

        #region Integrals
        public static Complex LineIntegralAlongCurve(this Func<Complex, Complex> function, Curve curve, int pointsCount = 10_000)
        {
            var sum = Complex.Zero;
            var points = curve.EvaluatePoints(pointsCount + 1);
            var result = curve.ApplyFunction(function).EvaluatePoints(pointsCount + 1);
            for (int point = 0; point < pointsCount; point++)
                sum += result[point] * (points[point + 1] - points[point]);
            return sum;
        }

        public static Task<Complex> LineIntegralAlongCurveAsync(this Func<Complex, Complex> function, Curve curve, int pointsCount = 10_000)
        {
            var job = new Task<Complex>(() =>
            {
                var sum = Complex.Zero;
                var points = curve.EvaluatePointsAsync(pointsCount + 1);
                points.Wait();
                var result = curve.ApplyFunction(function).EvaluatePointsAsync(pointsCount + 1);
                result.Wait();
                for (int point = 0; point < pointsCount; point++)
                    sum += result.Result[point] * (points.Result[point + 1] - points.Result[point]);
                return sum;
            });
            job.Start();
            return job;
        }
        #endregion
    }
}
