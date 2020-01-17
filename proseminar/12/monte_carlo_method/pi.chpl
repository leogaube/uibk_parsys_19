use Random;
use Time;

proc main(args: [] string){
	var count: sync int = 0; // number of points in the 1st quadrant of unit circle
	var N: int = 10000;
	if args.size >= 1 {
		N = args(1): int;
	}

	var timer = new Timer();
	timer.start();

	var randStream = new owned RandomStream(real, 0);
	forall i in 0..(N-1) {
		var x: real = randStream.getNext();
		var y: real = randStream.getNext();
		var z: real = x * x + y * y;
		if z <= 1 {
			count += 1;
		}
	}

	var pi: real = count / N:real * 4;

	timer.stop();
	writeln("calculated value of pi:",pi);
	writeln("The process took ",timer.elapsed(TimeUnits.seconds)," seconds to finish");

	return;
}
