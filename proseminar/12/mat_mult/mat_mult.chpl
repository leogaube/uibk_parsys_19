
use Time;
proc main(args: [] string) {
  var timer = new Timer();
  timer.start();
  
  if (args.size != 2) {
    writeln("Invalid number of arguments, please specify how many elements the arrays should have.");
    return;
  }

  var n: int;
  n = args(1): int;
  var seed = 17;
  var A: [1..n][1..n] int;
  var B: [1..n][1..n] int;
  var C: [1..n][1..n] int;

  // initialization
  for i in 1..n do {
    for j in 1..n do {
      A(i)(j) = j;
    }
    B(i)(i) = 1;
  }

  // Calcuation
  forall i in 1..n do {
    for j in 1..n do {
      for k in 1..n do {
        C(i)(j) += A(i)(k) * B(k)(j);
      }
    }
  }

  timer.stop();

  // Validation
  for i in 1..n do {
    for j in 1..n do {
      if (A(i)(j) != C(i)(j)) {
        writeln("Validation error, calcuation was not right!");
        return;
      }
    }
  }
  writeln("Validation succeeded!");

  writeln("Execution time is: ", timer.elapsed(TimeUnits.seconds));
  return;
}