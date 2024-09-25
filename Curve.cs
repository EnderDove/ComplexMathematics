namespace ComplexMathematics
{
    public class Curve(Func<double, Complex> parametricFunction, double startPoint, double endPoint)
    {
        public readonly Func<double, Complex> ParametricFunction = parametricFunction;
        public readonly double StartPoint = startPoint;
        public readonly double EndPoint = endPoint;

        public Complex[] EvaluatePoints(int pointCount)
        {
            var points = new Complex[pointCount + 1];
            double length = EndPoint - StartPoint;
            for (int i = 0; i < pointCount + 1; i++)
            {
                var t = StartPoint + (double)i / pointCount * length;
                points[i] = ParametricFunction(t);
            }
            return points;
        }
        public Task<Complex[]> EvaluatePointsAsync(int pointCount)
        {
            var job = new Task<Complex[]>(()=>
            {
                var points = new Complex[pointCount + 1];
                double length = EndPoint - StartPoint;
                Parallel.For(0, pointCount + 1, (int i) =>
                {
                    var t = StartPoint + (double)i / pointCount * length;
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

        public static Curve CreateFromPolarCoordinates(Func<double, Complex> polarFunction, double startPoint, double endPoint)
        {
            return new(ConvertToParametricFunction, startPoint, endPoint);

            Complex ConvertToParametricFunction(double phi)
            {
                var value = (polarFunction(phi));
                return new(value.Real * Math.Cos(phi), value.Imaginary * Math.Sin(phi));
            }
        }
    }
}
