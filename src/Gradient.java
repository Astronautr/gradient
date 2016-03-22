
import java.util.Arrays;
import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.interval.set.SetIntervalOps;
import net.java.jinterval.rational.ExtendedRational;
import net.java.jinterval.rational.Rational;

public class Gradient {

    /**
     * Gradient class for interval computing on Java. With this class you can
     * compute partial derivatives and the range of values for some functions.
     * Here was used the idea of automatic differentiation. JInterval library
     * was used for all interval arithmetic and functions.
     */
    private final Factory factory;
    private SetInterval X; //Interval value of expression.
    private final SetInterval dX[]; //Interval values of derivative.

    private static final SetInterval ZERO = SetIntervalOps.nums2(0, 0);
    private static final SetInterval ONE = SetIntervalOps.nums2(1, 1);

    /**
     * Owner of Gradients related to specific problem
     */
    public static class Factory {

        private final Gradient[] origin; //Represents the number of variables.
        private SetIntervalContext ic; //This field determine accuracy of computing.

        private Factory(SetIntervalContext ic, SetInterval x[]) {
            int dim = x.length;
            origin = new Gradient[dim];
            for (int i = 0; i < dim; i++) {
                Gradient xi = new Gradient(this);
                xi.X = x[i];
                if (xi.X == null) {
                    throw new NullPointerException();
                }
                xi.dX[i] = ONE;
                origin[i] = xi;
            }
            this.ic = ic;
        }

        /**
         * Return number of inputs in the problem.
         *
         * @return
         */
        public int getDim() {
            return origin.length;
        }

        /**
         * Return gradient for a given input variable
         *
         * @param i index of input variable starting from 0.
         * @return Gradient for an input variable
         */
        public Gradient getInp(int i) {
            return origin[i];
        }

        /**
         * Return array of Gradients for all input variables
         *
         * @return Gradients
         */
        public Gradient[] getInps() {
            return origin.clone();
        }

        /**
         * Return Gradient for a constant
         *
         * @param number value of a constant
         * @return Gradient for a constant
         */
        public Gradient num(double number) {
            Gradient result = new Gradient(this);
            result.X = ic.numsToInterval(number, number);
            return result;
        }

        /**
         * Return Gradient for a constant
         *
         * @param number value of a constant
         * @return Gradient for a constant
         */
        public Gradient num(Rational number) {
            Gradient result = new Gradient(this);
            result.X = ic.numsToInterval(number, number);
            return result;
        }

        /**
         * Return Gradient for an unknown constant
         *
         * @param lower lower bound of constant
         * @param upper upper bound of constant
         * @return Gradient for an unknown constant
         */
        public Gradient nums(double lower, double upper) {
            Gradient result = new Gradient(this);
            result.X = ic.numsToInterval(lower, upper);
            return result;
        }

        /**
         * Return Gradient for an unknown constant
         *
         * @param lower lower bound of constant
         * @param upper upper bound of constant
         * @return Gradient for an unknown constant
         */
        public Gradient nums(Rational lower, Rational upper) {
            Gradient result = new Gradient(this);
            result.X = ic.numsToInterval(lower, upper);
            return result;
        }

        /**
         * Return Gradient for an unknown constant
         *
         * @param lower lower bound of constant
         * @param upper upper bound of constant
         * @return Gradiend for an unknown constant
         */
        public Gradient nums(ExtendedRational lower, ExtendedRational upper) {
            Gradient result = new Gradient(this);
            result.X = ic.numsToInterval(lower, upper);
            return result;
        }

        private void check(Gradient that) {
            if (that.factory != this) {
                throw new IllegalArgumentException();
            }
        }
    }

    /**
     * Create gradient Factory for a problem. Problem dimension is specified by
     * a box.
     *
     * @param ic interval context for internal computatios
     * @param x input box.
     * @return Gradient Factory for a new problem
     */
    public static Factory createFactory(SetIntervalContext ic, SetInterval... x) {
        return new Factory(ic, x);
    }

    private Gradient(Factory factory) {
        this.factory = factory;
        X = null;
        dX = new SetInterval[factory.getDim()];
        Arrays.fill(dX, ZERO);
    }

    private void copyGrads(Gradient that) {
        System.arraycopy(this.dX, 0, that.dX, 0, dX.length);
    }

    /**
     * Factory to which this Gradient belongs
     *
     * @return Factory of this Gradient
     */
    public Factory getFactory() {
        return factory;
    }

    /**
     * Demension of the problem
     *
     * @return dimension of the problem
     */
    public int getDim() {
        return dX.length;
    }

    public Gradient neg() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.neg(X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.neg(dX[i]);
        }
        return result;
    }

    public Gradient intersection(Gradient Y) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.intersection(this.X, Y.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.intersection(this.dX[i], Y.dX[i]);
        }
        return result;
    }

    public Gradient intersectionX(Gradient Y) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.intersection(this.X, Y.X);
        result.copyGrads(this);
        return result;
    }

    public Gradient intersectionX(SetInterval Y) {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.intersection(this.X, Y);
        result.copyGrads(this);
        return result;
    }

    public Gradient intersectionDX(Gradient Y) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = this.X;
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.intersection(this.dX[i], Y.dX[i]);
        }
        return result;
    }

    public Gradient intersectionDX(SetInterval[] Y) {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = this.X;
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.intersection(this.dX[i], Y[i]);
        }
        return result;
    }

    public Gradient intersectionDXInd(Gradient Y, int ind) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = this.X;
        result.copyGrads(this);
        result.dX[ind] = ic.intersection(this.dX[ind], Y.dX[ind]);
        return result;
    }

    public Gradient intersectionDXInd(SetInterval Y, int ind) {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = this.X;
        result.copyGrads(this);
        result.dX[ind] = ic.intersection(this.dX[ind], Y);
        return result;
    }

    public Gradient intersectionDXInd(SetInterval[] Y, int ind) {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = this.X;
        result.copyGrads(this);
        result.dX[ind] = ic.intersection(this.dX[ind], Y[ind]);
        return result;
    }

    public Gradient hull(Gradient Y) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.convexHull(this.X, Y.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.convexHull(this.dX[i], Y.dX[i]);
        }
        return result;
    }

    public Gradient add(Gradient Y) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.add(this.X, Y.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.add(this.dX[i], Y.dX[i]);
        }
        return result;
    }

    public Gradient sub(Gradient Y) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.sub(this.X, Y.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.sub(this.dX[i], Y.dX[i]);
        }
        return result;
    }

    public Gradient mul(Gradient Y) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.mul(this.X, Y.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.add(ic.mul(this.X, Y.dX[i]), ic.mul(this.dX[i], Y.X));
        }
        return result;
    }

    public Gradient div(Gradient Y) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.div(this.X, Y.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.div(ic.sub(ic.mul(this.dX[i], Y.X), ic.mul(this.X, Y.dX[i])), ic.sqr(Y.X));
        }
        return result;
    }

    public Gradient pow(Gradient Y) {
        factory.check(Y);
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.pow(this.X, Y.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.mul(ic.mul(Y.X, ic.pow(X, ic.sub(Y.X, ic.numsToInterval(1, 1)))), this.dX[i]);
        }
        return result;
    }

    public Gradient pown(int n) {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.pown(this.X, n);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(n, n), ic.pown(X, n - 1)), this.dX[i]);
        }
        return result;
    }

    public Gradient pown(long n) {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.pown(this.X, n);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(n, n), ic.pown(X, n - 1)), this.dX[i]);
        }
        return result;
    }

    public Gradient sqr() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.sqr(this.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(2, 2), X), this.dX[i]);
        }
        return result;
    }

    public Gradient sqrt() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.sqrt(this.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.div(dX[i], ic.mul(ic.numsToInterval(2, 2), ic.sqrt(X)));
        }
        return result;
    }

    public Gradient sin() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.sin(this.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.mul(ic.cos(X), dX[i]);
        }
        return result;
    }

    public Gradient cos() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.cos(this.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.neg(ic.mul(ic.sin(X), dX[i]));
        }
        return result;
    }

    public Gradient tan() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.tan(this.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.mul(ic.div(ic.numsToInterval(1, 1), ic.sqr(ic.cos(X))), dX[i]);
        }
        return result;
    }

    public Gradient asin() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.asin(this.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.mul(ic.div(ic.numsToInterval(1, 1), ic.sqrt(ic.sub(ic.numsToInterval(1, 1), ic.sqr(X)))), dX[i]);
        }
        return result;
    }

    public Gradient acos() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.acos(this.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.neg(ic.mul(ic.div(ic.numsToInterval(1, 1), ic.sqrt(ic.sub(ic.numsToInterval(1, 1), ic.sqr(X)))), dX[i]));
        }
        return result;
    }

    public Gradient atan() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.atan(this.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.mul(ic.div(ic.numsToInterval(1, 1), ic.add(ic.numsToInterval(1, 1), ic.sqr(X))), dX[i]);
        }
        return result;
    }

    public Gradient exp() {
        SetIntervalContext ic = factory.ic;
        Gradient result = new Gradient(factory);
        result.X = ic.exp(this.X);
        for (int i = 0; i < getDim(); i++) {
            result.dX[i] = ic.mul(ic.exp(X), dX[i]);
        }
        return result;
    }

    public SetInterval getX() {
        return X;
    }

    public SetInterval getDX(int i) {
        return dX[i];
    }

    public SetInterval[] getDX() {
        return dX.clone();
    }

    public void show() {
        System.out.println("[" + this.X.doubleInf() + ", " + this.X.doubleSup() + "]");
        System.out.print("(");
        for (int i = 0; i < getDim(); i++) {
            System.out.print(" [" + this.dX[i].doubleInf() + ", " + this.dX[i].doubleSup() + "]");
        }
        System.out.println(" )");
    }
}
