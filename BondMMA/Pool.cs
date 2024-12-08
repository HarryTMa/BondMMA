namespace BondMMA
{
    public class Pool(double y0, double r0 = 0.05, double k = 0.02)
    {
        /// <summary>
        /// The current value of the bonds in the pool. The initial value is y0 according to the derivation
        /// </summary>
        public double X { get; private set; } = y0;
        /// <summary>
        /// The cash balance in the pool
        /// </summary>
        public double y { get; private set; } = y0;
        /// <summary>
        /// kappa, the parameter
        /// </summary>
        public double k { get; private set; } = k;
        /// <summary>
        /// r*, the anchor rate
        /// </summary>
        public double R { get; set; } = r0;
        private readonly double y0 = y0;
        double loans = 0;
        public double NetEquity => y - y0 + loans;
        public double? TransactByY(double dy, double t, double? R = null)
        {
            if (dy > 0 && NetEquity < y0 * -0.01)
            {
                return null;
            }
            double a = 1 / (1 + k * t);
            double dx = Math.Exp((R ?? this.R) * t) * y * (Math.Pow(X / y + 1 - Math.Pow(dy / y + 1, a), 1 / a) - Math.Pow(X / y, 1 / a));
            double newY = y + dy;
            X = Math.Pow(y, a - 1) * (X + y) / Math.Pow(newY, a - 1) - newY;
            y = newY;
            loans -= dy;
            return dx;
        }
        public double TransactByX(double dx, double t, double? R = null)
        {
            double a = 1 / (1 + k * t);
            double dy = y * (Math.Pow(X / y + 1 - Math.Pow(Math.Pow(X / y, 1 / a) + Math.Exp(-(R ?? this.R) * t) * dx / y, a), 1 / a) - 1);
            double newY = y + dy;
            X = Math.Pow(y, a - 1) * (X + y) / Math.Pow(newY, a - 1) - newY;
            y = newY;
            loans -= dy;
            return dy;
        }
        /// <summary>
        /// The margin rate
        /// </summary>
        /// <param name="R">The anchor rate</param>
        /// <returns>The margin rate</returns>
        public double r(double? R = null)
        {
            return k * Math.Log(X / y) + (R ?? this.R);
        }
        public void Proceed(double dt)
        {
            loans *= Math.Exp(dt * r());
        }
    }
}
