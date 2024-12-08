namespace BondMMA
{
    internal static class Program
    {
        public static double Normal(this Random rnd, double mean = 0, double std = 1)
        {
            double u1 = 1.0 - rnd.NextDouble();
            double u2 = 1.0 - rnd.NextDouble();
            double z = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
            return mean + std * z;
        }

        static void Main(string[] args)
        {
            Random rnd = new(Seed: int.Parse(args[0]));
            double r = 0.05; // r_0 = 0.05
            Pool pool = new(1000, r);
            int N = 100000;
            int M = 1000;
            double k = 0.4;
            double theta = 0.05;
            double sigma = .2;
            double dt = 1.0 / N;
            double[] market_r = new double[N];
            double[] pool_r = new double[N];
            double[] pool_r_i = new double[M];
            double[] diff = new double[N];
            double[] pool_r_std = new double[N];
            double[] equity = new double[N];
            double[] reject = new double[N];
            List<double>[] loans = new List<double>[N];
            double cash = 0;
            for (int i = 0; i < N; i++)
            {
                double T = i / (double)N;
                r += k * (theta - r) * dt + sigma * Math.Sqrt(r * dt) * rnd.Normal();
                market_r[i] = r;
                // Pay off the loans
                double P;
                double p = 0;
                if (loans[i] != null)
                {
                    P = loans[i].Count / (double)M;
                }
                else
                {
                    P = 0;
                }
                int l = -1;
                for (int j = 0; j < M; j++)
                {
                    p += P;
                    while (p > 1)
                    {
                        p -= 1;
                        l += 1;
                        double dx_ = loans[i][l];
                        double dy_ = pool.TransactByX(-dx_, 0);
                        cash -= dy_;
                    }
                    double r_ = pool.r();
                    pool_r_i[j] = r_;
                    double dy;
                    if (r_ < r) // Should Borrow
                    {
                        dy = -Math.Abs(rnd.Normal(0.72)); // mean is approximately 1
                    }
                    else
                    {
                        dy = Math.Abs(rnd.Normal(0.72));
                    }
                    double t = Math.Abs(rnd.Normal(1 - T, 1 - T));
                    int ticks = (int)(t * N) + 1;
                    t = ticks * dt;
                    ticks += i;
                    double? dx = pool.TransactByY(dy, t);
                    if (dx == null)
                    {
                        reject[i] += 1.0 / M;
                        continue;
                    }
                    cash -= dy;
                    if (ticks < N) // The loan will be paid off
                    {
                        if (loans[ticks] == null)
                        {
                            loans[ticks] = [];
                        }
                        loans[ticks].Add(dx.Value);
                    }
                }
                pool.Proceed(dt);
                pool_r[i] = pool_r_i.Sum() / M;
                pool_r_std[i] = Math.Sqrt(pool_r_i.Select(x => (x - pool_r[i]) * (x - pool_r[i])).Sum() / M);
                equity[i] = pool.NetEquity;
                diff[i] = pool_r[i] - market_r[i];
                pool.R = r;
            }
            ScottPlot.Plot plot = new();
            double[] indices = Enumerable.Range(0, N).Select(i => i / (double)N).ToArray();
            plot.Add.Scatter(indices, pool_r, ScottPlot.Color.FromColor(System.Drawing.Color.Orange));
            plot.Add.Scatter(indices, market_r, ScottPlot.Color.FromColor(System.Drawing.Color.Blue));
            plot.SavePng("r.png", 1000, 1000);
            plot = new();
            plot.Add.Scatter(indices, diff);
            plot.SavePng("diff.png", 1000, 1000);
            plot = new();
            plot.Add.Scatter(indices, equity);
            plot.SavePng("equity.png", 1000, 1000);
            plot = new();
            plot.Add.Scatter(indices, pool_r_std);
            plot.SavePng("std.png", 1000, 1000);
            plot = new();
            plot.Add.Scatter(indices, reject);
            plot.SavePng("reject.png", 1000, 1000);
            Console.WriteLine("Cash remains: {0}", cash);
        }
    }
}
