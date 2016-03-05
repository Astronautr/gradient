import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.rational.BinaryValueSet;


public class Main {
    //Example of using Gradient with JInterval
    public static void main(String[] args) {
        SetIntervalContext ic = SetIntervalContexts.getFast();
        System.out.println("First example. getFast context was used.");
        SetInterval[] box = {ic.numsToInterval(1,2), ic.numsToInterval(2,3), ic.numsToInterval(0.2,1.18)};
        Gradient[] origin = Gradient.init(box, ic);
        System.out.println("And the first expression is: f(x[]) = x[0]^2 + (x[1] + 1)^2 - 2*x[0]*(x[1] + 1) * sin(x[2])");
        (origin[0].sqr().add((origin[1].add(Gradient.num(1))).sqr()).
                sub(origin[0].mul(Gradient.num(2)).mul(origin[1].add(Gradient.num(1))).
                        mul(origin[2].sin()))).show();

        System.out.println("\nNow we change the accuracy and dimension of example");
        ic = SetIntervalContexts.getInfSup(BinaryValueSet.BINARY64);
        Gradient x = new Gradient(ic.numsToInterval(0.999,1.001), ic);
        System.out.println("f(x) = sin(x) * (4*cos(x) - 2)^2");
        x.sin().mul(
                (Gradient.num(4).mul(x.cos()).sub(Gradient.num(2))).sqr()
        ).show();

        System.out.println("\nNow we take more difficult expression");
        System.out.println("Here we have two functions:" +
                "\ng(x[]) =sqrt(x[0]^2 + x[1]^4) * exp(-x[2])" +
                "\nf(x[],y) = y / (x[2]^3 + atan(x[0] * x[1]))" +
                "\n Find f(x[], g(x[]))");
        box = new SetInterval[]{ic.numsToInterval(1,2), ic.numsToInterval(4,7), ic.numsToInterval(1,5)};
        origin = Gradient.init(box,ic);
        System.out.println("\nComputing g(x[]).");
        Gradient g = (origin[0].sqr().add(origin[1].pown(4))).sqrt().mul((origin[2].neg()).exp());
        g.show();
        System.out.println("\nComputing f(x[], g(x[])).");
        Gradient f = g.div(origin[2].pown(3).add((origin[0].mul(origin[1])).atan()));
        f.show();
    }
}
