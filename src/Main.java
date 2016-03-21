
import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.rational.BinaryValueSet;

public class Main {

    //Example of using Gradient with JInterval
    public static void main(String[] args) {
        SetIntervalContext ic = SetIntervalContexts.getPlain();
        System.out.println("First example. getFast context was used.");
        SetInterval[] box = {ic.numsToInterval(1, 2), ic.numsToInterval(2, 3), ic.numsToInterval(0.2, 1.18)};
        Gradient.Factory fa = Gradient.createFactory(ic, box);
        Gradient x0 = fa.getInp(0), x1 = fa.getInp(1), x2 = fa.getInp(2);
        System.out.println("And the first expression is: f(x[]) = x[0]^2 + (x[1] + 1)^2 - 2*x[0]*(x[1] + 1) * sin(x[2])");
        (x0.sqr().add((x1.add(fa.num(1))).sqr()).
                sub(x0.mul(fa.num(2)).mul(x1.add(fa.num(1))).
                        mul(x2.sin()))).show();

        System.out.println("\nNow we change the accuracy and dimension of example");
        ic = SetIntervalContexts.getInfSup(BinaryValueSet.BINARY64);
        Gradient.Factory fb = Gradient.createFactory(ic, ic.numsToInterval(0.999, 1.001));
        Gradient x = fb.getInp(0);
        System.out.println("f(x) = sin(x) * (4*cos(x) - 2)^2");
        x.sin().mul((fb.num(4).mul(x.cos()).sub(fb.num(2))).sqr()
        ).show();

        System.out.println("\nNow we take more difficult expression");
        System.out.println("Here we have two functions:"
                + "\ng(x[]) =sqrt(x[0]^2 + x[1]^4) * exp(-x[2])"
                + "\nf(x[],y) = y / (x[2]^3 + atan(x[0] * x[1]))"
                + "\n Find f(x[], g(x[]))");
        ic = SetIntervalContexts.getAccur64();
        box = new SetInterval[]{ic.numsToInterval(1, 2), ic.numsToInterval(4, 7), ic.numsToInterval(1, 5)};
        Gradient.Factory fc = Gradient.createFactory(ic, box);
        x0 = fc.getInp(0);
        x1 = fc.getInp(1);
        x2 = fc.getInp(2);
        System.out.println("\nComputing g(x[]).");
        Gradient g = (x0.sqr().add(x1.pown(4))).sqrt().mul((x2.neg()).exp());
        g.show();
        System.out.println("\nComputing f(x[], g(x[])).");
        Gradient f = g.div(x2.pown(3).add((x0.mul(x1)).atan()));
        f.show();
    }
}
