import net.java.jinterval.interval.set.SetInterval;
import net.java.jinterval.interval.set.SetIntervalContext;
import net.java.jinterval.interval.set.SetIntervalContexts;
import net.java.jinterval.rational.BinaryValueSet;

/**
 * Created by astronaut on 3/2/16.
 */
public class Main {
    //Example of using Gradient with JInterval
    public static void main(String[] args) {
        SetIntervalContext ic = SetIntervalContexts.getFast();
        SetInterval[] bar = {ic.numsToInterval(1,2), ic.numsToInterval(2,3), ic.numsToInterval(0,1.18)};
        Gradient[] origin = Gradient.init(bar, ic);
        (origin[0].mul(origin[0]).add((origin[1].add(Gradient.num(1))).sqr()).
                sub(origin[0].mul(Gradient.num(2)).mul(origin[1].add(Gradient.num(1))).
                        mul(origin[2].sin()))).show();
        origin[0].hull(origin[1]).show();

    }
}
