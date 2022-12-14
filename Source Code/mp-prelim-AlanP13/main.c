#include "hello.h"

/******************************************************************************
 * 
 * Prelim lab header
 * 
 * Author: Alan Palayil
 * Email:  apalayil@hawk.iit.edu
 * AID:    A20447935
 * Date:   8/30/21
 * 
 * By signing above, I pledge on my honor that I neither gave nor received any
 * unauthorized assistance on the code contained in this repository.
 * 
 *****************************************************************************/

int main(int argc, char *argv[]) {
  if (argc > 1) {
    say_hello_to(argv[1]);
  } else {
    say_hello_to("world");
  }
  return 0;
}
