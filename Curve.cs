namespace ComplexMathematics
{
    /// <summary>
    /// A Class represents Curve on the complex plane given by <paramref name="parametricFunction"/> with <paramref name="startPoint"/> and <paramref name="endPoint"/>.
    /// </summary>
    public class Curve(Func<double, Complex> parametricFunction, double startPoint, double endPoint)
    {
        public readonly Func<double, Complex> ParametricFunction = parametricFunction;
        public readonly double StartPoint = startPoint;
        public readonly double EndPoint = endPoint;

        public static Curve CreateFromPolarCoordinates(Func<double, double> polarFunction, double startPoint, double endPoint)
        {
            return new(ConvertToParametricFunction, startPoint, endPoint);

            Complex ConvertToParametricFunction(double phi)
            {
                var value = polarFunction(phi);
                return new(value * Math.Cos(phi), value * Math.Sin(phi));
            }
        }

        public Complex[] EvaluatePoints(int pointCount)
        {
            var points = new Complex[pointCount];
            double length = EndPoint - StartPoint;
            for (int i = 0; i < pointCount; i++)
            {
                var t = StartPoint + i / (pointCount - 1.0) * length;
                points[i] = ParametricFunction(t);
            }
            return points;
        }
        public Task<Complex[]> EvaluatePointsAsync(int pointCount)
        {
            var job = new Task<Complex[]>(() =>
            {
                var points = new Complex[pointCount];
                double length = EndPoint - StartPoint;
                Parallel.For(0, pointCount, (int i) =>
                {
                    var t = StartPoint + i / (pointCount - 1.0) * length;
                    points[i] = ParametricFunction(t);
                });
                return points;
            });
            job.Start();
            return job;
        }

        public Curve ApplyFunction(Func<Complex, Complex> function)
        {
            return new((double value) => function(ParametricFunction(value)), StartPoint, EndPoint);
        }
    }
}
