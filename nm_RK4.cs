using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MetaPop_SIR_Model
{
    internal class nm_RK4
    {
        double m_S1, m_I1, m_R1, P1, m_S2, m_I2, m_R2, P2;
        double m_BETA1, m_GAMMA1, m_BETA2, m_GAMMA2, m_PHI1, m_THETA1, m_PHI2, m_THETA2;
        int NLast;
        double h;

        public nm_RK4(double P10, double P20,
            double S10, double S20,
            double I10, double I20,
            double R10, double R20,
            double BETA1, double BETA2, 
            double GAMMA1, double GAMMA2,
            double PHI1, double PHI2,
            double THETA1, double THETA2,
            int DT, 
            double SS)
        {
            P1 = P10; P2 = P20;
            m_S1 = S10 / P10; m_S2 = S20 / P20;
            m_I1 = I10 / P10; m_I2 = I20 / P20;
            m_R1 = R10 / P10; m_R2 = R20 / P20;
            m_BETA1 = BETA1; m_BETA2 = BETA2;
            m_GAMMA1 = GAMMA1; m_GAMMA2 = GAMMA2;
            m_PHI1 = PHI1; m_PHI2 = PHI2;
            m_THETA1 = THETA1; m_THETA2 = THETA2;
            NLast = DT;
            h = SS;
        }

        public void fluSim(double[,] suS, double[,] inF, double[,] reM,
                           double[] PHI, double[] THETA)
        {
            double S1, S2, S3, S4, SA, I1, I2, I3, I4, IA, R1, R2, R3, R4, RA;
            double A, B;
            double t = 0.0;
            suS[0, 0] = m_S1; suS[0, 1] = m_S2;
            inF[0, 0] = m_I1; inF[0, 1] = m_I2;
            reM[0, 0] = m_R1; reM[0, 1] = m_R2;
            PHI[0] = m_PHI1; PHI[1] = m_PHI2;
            THETA[0] = m_THETA1; THETA[1] = m_THETA2;
            double[] TRANS = new double[2];

            for (int i = 1; i <= (NLast / h); i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    TRANS[j] = PHI[j] * Math.Exp(-THETA[j] * t);
                    if (j == 0)
                    {
                        A = TRANS[0] * inF[i - 1, 0];
                        B = TRANS[1] * inF[i - 1, 1];
                        //Correction
                        S1 = h * f1(j, t, m_S1, m_I1, m_R1);
                        I1 = h * f2(j, t, m_S1, m_I1, m_R1, A, B);
                        R1 = h * f3(j, t, m_S1, m_I1, m_R1);

                        S2 = h * f1(j, t + h / 2, m_S1 + S1 / 2, m_I1 + I1 / 2, m_R1 + R1 / 2);
                        I2 = h * f2(j, t + h / 2, m_S1 + S1 / 2, m_I1 + I1 / 2, m_R1 + R1 / 2, A, B);
                        R2 = h * f3(j, t + h / 2, m_S1 + S1 / 2, m_I1 + I1 / 2, m_R1 + R1 / 2);

                        S3 = h * f1(j, t + h / 2, m_S1 + S2 / 2, m_I1 + I2 / 2, m_R1 + R2 / 2);
                        I3 = h * f2(j, t + h / 2, m_S1 + S2 / 2, m_I1 + I2 / 2, m_R1 + R2 / 2, A, B);
                        R3 = h * f3(j, t + h / 2, m_S1 + S2 / 2, m_I1 + I2 / 2, m_R1 + R2 / 2);

                        S4 = h * f1(j, t + h, m_S1 + S3, m_I1 + I3, m_R1 + R3);
                        I4 = h * f2(j, t + h, m_S1 + S3, m_I1 + I3, m_R1 + R3, A, B);
                        R4 = h * f3(j, t + h, m_S1 + S3, m_I1 + I3, m_R1 + R3);

                        SA = 1.0 / 6.0 * (S1 + 2 * S2 + 2 * S3 + S4);
                        IA = 1.0 / 6.0 * (I1 + 2 * I2 + 2 * I3 + I4);
                        RA = 1.0 / 6.0 * (R1 + 2 * R2 + 2 * R3 + R4);

                        m_S1 = m_S1 + SA;
                        m_I1 = m_I1 + IA;
                        m_R1 = m_R1 + RA;

                        suS[i, 0] = m_S1;
                        inF[i, 0] = m_I1;
                        reM[i, 0] = m_R1;
                    }
                    else
                    {
                        A = TRANS[0] * inF[i, 0];
                        B = TRANS[1] * inF[i - 1, 1];
                        //Correction
                        S1 = h * f1(j, t, m_S2, m_I2, m_R2);
                        I1 = h * f2(j, t, m_S2, m_I2, m_R2, A, B);
                        R1 = h * f3(j, t, m_S2, m_I2, m_R2);

                        S2 = h * f1(j, t + h / 2, m_S2 + S1 / 2, m_I2 + I1 / 2, m_R2 + R1 / 2);
                        I2 = h * f2(j, t + h / 2, m_S2 + S1 / 2, m_I2 + I1 / 2, m_R2 + R1 / 2, A, B);
                        R2 = h * f3(j, t + h / 2, m_S2 + S1 / 2, m_I2 + I1 / 2, m_R2 + R1 / 2);

                        S3 = h * f1(j, t + h / 2, m_S2 + S2 / 2, m_I2 + I2 / 2, m_R2 + R2 / 2);
                        I3 = h * f2(j, t + h / 2, m_S2 + S2 / 2, m_I2 + I2 / 2, m_R2 + R2 / 2, A, B);
                        R3 = h * f3(j, t + h / 2, m_S2 + S2 / 2, m_I2 + I2 / 2, m_R2 + R2 / 2);

                        S4 = h * f1(j, t + h, m_S2 + S3, m_I2 + I3, m_R2 + R3);
                        I4 = h * f2(j, t + h, m_S2 + S3, m_I2 + I3, m_R2 + R3, A, B);
                        R4 = h * f3(j, t + h, m_S2 + S3, m_I2 + I3, m_R2 + R3);

                        SA = 1.0 / 6.0 * (S1 + 2 * S2 + 2 * S3 + S4);
                        IA = 1.0 / 6.0 * (I1 + 2 * I2 + 2 * I3 + I4);
                        RA = 1.0 / 6.0 * (R1 + 2 * R2 + 2 * R3 + R4);

                        m_S2 = m_S2 + SA;
                        m_I2 = m_I2 + IA;
                        m_R2 = m_R2 + RA;

                        suS[i, 1] = m_S2;
                        inF[i, 1] = m_I2;
                        reM[i, 1] = m_R2;
                    }
                }
                t = t + h;  
            }
        }

        private double f1(int j, double t, double S, double D, double R)
        {
            if (j == 0)
            {
                return -m_BETA1 * P1 * S * D;
            }
            else
            {
                return -m_BETA2 * P2 * S * D;
            }

        }

        private double f2(int j, double t, double S, double D, double R, double A, double B)
        {
            if (j == 0)
            {
                return m_BETA1 * P1 * S * D - m_GAMMA1 * D - A + B;
            }
            else
            {
                return m_BETA2 * P2 * S * D - m_GAMMA2 * D + A - B;
            }
        }

        private double f3(int j, double t, double S, double D, double R)
        {
            if (j == 0)
            {
                return m_GAMMA1 * D;
            }
            else
            {
                return m_GAMMA2 * D;
            }
        }
    }
}
