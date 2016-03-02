import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.rational.BinaryValueSet;

/**
 * Created by astronaut on 2/22/16.
 */
public class Gradient {
    private static SetIntervalContext ic;
    private static int dim;
    private SetInterval X;
    private SetInterval dX[];

    public Gradient() {
        //ic = SetIntervalContexts.getInfSup(BinaryValueSet.BINARY64);
        //ic = SetIntervalContexts.getFast();
        //dim = 0;
        X = null;
        dX = null;
    }

    public Gradient(SetInterval bar, SetIntervalContext ic) {
        dim = 1;
        dX = new SetInterval[1];
        dX[0] = ic.numsToInterval(1,1);
        X = bar;
        this.ic = ic;
    }

    public static Gradient[] init(SetInterval[] bar, SetIntervalContext ic) {
        dim = bar.length;
        Gradient.ic = ic;
        Gradient[] result = new Gradient[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = new Gradient();
            result[i].dX = new SetInterval[dim];
            result[i].X = bar[i];
            for (int j = 0; j < dim; j++) { //Здесь могло быть dim сравнений и dim присваиваний
                result[i].dX[j] = ic.numsToInterval(0,0);
            }
            result[i].dX[i] = ic.numsToInterval(1,1);
        }
        return result;
    }

    public static Gradient num(double number) {
        Gradient result = new Gradient();
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.numsToInterval(0,0);
        }
        result.X = ic.numsToInterval(number,number);
        return result;
    }

    public Gradient neg() {
        Gradient result = new Gradient();
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(dX[i]);
        }
        result.X = ic.neg(X);
        return result;
    }

    public Gradient intersection(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.intersection(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.intersection(this.dX[i],Y.dX[i]);
        }
        return result;
    }

    public Gradient intersectionX(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.intersection(this.X, Y.X);
        result.dX = new SetInterval[dim];
        result.dX = this.dX.clone();
        return result;
    }

    public Gradient intersectionX(SetInterval Y) {
        Gradient result = new Gradient();
        result.X = ic.intersection(this.X, Y);
        result.dX = new SetInterval[dim];
        result.dX = this.dX.clone();
        return result;
    }

    public Gradient add(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.add(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.add(this.dX[i],Y.dX[i]);
        }
        return result;
    }

    public Gradient sub(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.sub(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.sub(this.dX[i],Y.dX[i]);
        }
        return result;
    }

    public Gradient mul(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.mul(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.add(ic.mul(this.X, Y.dX[i]),ic.mul(this.dX[i],Y.X));
        }
        return result;
    }

    public Gradient div(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.div(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(ic.sub(ic.mul(this.dX[i], Y.X),ic.mul(this.X,Y.dX[i])),ic.sqr(Y.X));
        }
        return result;
    }

    public Gradient pow(Gradient Y) {
        Gradient result = new Gradient();
        result.X = ic.pow(this.X, Y.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(Y.X, ic.pow(X,ic.sub(Y.X, ic.numsToInterval(1,1)))),this.dX[i]);
        }
        return result;
    }

    public Gradient sqr() {
        Gradient result = new Gradient();
        result.X = ic.sqr(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.mul(ic.numsToInterval(2,2),X),this.dX[i]);
        }
        return result;
    }

    public Gradient sqrt() {
        Gradient result = new Gradient();
        result.X = ic.sqrt(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.div(dX[i], ic.mul(ic.numsToInterval(2,2),ic.sqrt(X)));
        }
        return result;
    }

    public Gradient sin() {
        Gradient result = new Gradient();
        result.X = ic.sin(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.cos(X),dX[i]);
        }
        return result;
    }

    public Gradient cos() {
        Gradient result = new Gradient();
        result.X = ic.cos(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.mul(ic.sin(X),dX[i]));
        }
        return result;
    }

    public Gradient tan() {
        Gradient result = new Gradient();
        result.X = ic.tan(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.div(ic.numsToInterval(1,1),ic.sqr(ic.cos(X))),dX[i]);
        }
        return result;
    }

    public Gradient asin() {
        Gradient result = new Gradient();
        result.X = ic.asin(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.div(ic.numsToInterval(1,1),ic.sqrt(ic.sub(ic.numsToInterval(1,1),ic.sqr(X)))),dX[i]);
        }
        return result;
    }

    public Gradient acos() {
        Gradient result = new Gradient();
        result.X = ic.acos(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.neg(ic.mul(ic.div(ic.numsToInterval(1,1),ic.sqrt(ic.sub(ic.numsToInterval(1,1),ic.sqr(X)))),dX[i]));
        }
        return result;
    }

    public Gradient atan() {
        Gradient result = new Gradient();
        result.X = ic.atan(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.div(ic.numsToInterval(1,1),ic.add(ic.numsToInterval(1,1),ic.sqr(X))),dX[i]);
        }
        return result;
    }

    public Gradient exp() {
        Gradient result = new Gradient();
        result.X = ic.exp(this.X);
        result.dX = new SetInterval[dim];
        for (int i = 0; i < dim; i++) {
            result.dX[i] = ic.mul(ic.exp(X),dX[i]);
        }
        return result;
    }

    public SetInterval getX() {
        return X;
    }

    public SetInterval[] getDX() {
        return dX;
    }

    public static double getDim() {return dim;}

    public void show() {
        System.out.println("[" + this.X.doubleInf() + ", " + this.X.doubleSup() + "]");
        System.out.print("(");
        for (int i = 0; i < dim; i++) {
            System.out.print(" [" + this.dX[i].doubleInf() + ", " + this.dX[i].doubleSup() + "]");
        }
        System.out.println(" )");
    }

}