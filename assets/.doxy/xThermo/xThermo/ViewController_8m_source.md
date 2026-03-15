

# File ViewController.m

[**File List**](files.md) **>** [**app\_ios**](dir_e936cb0e42e32aaf0ec84f84716b9a91.md) **>** [**ViewController.m**](ViewController_8m.md)

[Go to the documentation of this file](ViewController_8m.md)


```Objective-C
#import "ViewController.h"
#import "CppInterface.h"

@interface ViewController ()
{
    CppInterface* i;
}
@end

@implementation ViewController

- (void)viewDidLoad {
    [super viewDidLoad];
    i = [[CppInterface alloc]init];
}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

@end
```


