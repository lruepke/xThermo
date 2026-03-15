

# Group BITMASK\_phi



[**Modules**](modules.md) **>** [**BITMASK\_phi**](group__BITMASK__phi.md)





































































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**Update\_phi\_all**](group__BITMASK__phi.md#define-update_phi_all)  `[**Update\_phi\_d**](group__BITMASK__phi.md#define-update_phi_d)\|[**Update\_phi\_dd**](group__BITMASK__phi.md#define-update_phi_dd)\|[**Update\_phi\_t**](group__BITMASK__phi.md#define-update_phi_t)\|[**Update\_phi\_tt**](group__BITMASK__phi.md#define-update_phi_tt)\|[**Update\_phi\_dt**](group__BITMASK__phi.md#define-update_phi_dt)`<br> |
| define  | [**Update\_phi\_d**](group__BITMASK__phi.md#define-update_phi_d)  `2`<br> |
| define  | [**Update\_phi\_dd**](group__BITMASK__phi.md#define-update_phi_dd)  `4`<br> |
| define  | [**Update\_phi\_dt**](group__BITMASK__phi.md#define-update_phi_dt)  `32`<br> |
| define  | [**Update\_phi\_t**](group__BITMASK__phi.md#define-update_phi_t)  `8`<br> |
| define  | [**Update\_phi\_tt**](group__BITMASK__phi.md#define-update_phi_tt)  `16`<br> |

## Macro Definition Documentation





### define Update\_phi\_all 

```
#define Update_phi_all `Update_phi_d | Update_phi_dd | Update_phi_t | Update_phi_tt | Update_phi_dt`
```



Update all 


        

<hr>



### define Update\_phi\_d 

```
#define Update_phi_d `2`
```



\(2^1\): \(\left(\frac{\partial \phi}{\partial \delta} \right)_{\tau}\) 


        

<hr>



### define Update\_phi\_dd 

```
#define Update_phi_dd `4`
```



\(2^2\): \(\left(\frac{\partial^2 \phi}{\partial \delta^2} \right)_{\tau}\) 


        

<hr>



### define Update\_phi\_dt 

```
#define Update_phi_dt `32`
```



\(2^5\): \(\left(\frac{\partial^2 \phi}{\partial \delta \partial \tau} \right)\) 


        

<hr>



### define Update\_phi\_t 

```
#define Update_phi_t `8`
```



\(2^3\): \(\left(\frac{\partial \phi}{\partial \tau} \right)_{\delta}\) 


        

<hr>



### define Update\_phi\_tt 

```
#define Update_phi_tt `16`
```



\(2^4\): \(\left(\frac{\partial^2 \phi}{\partial \tau^2} \right)_{\delta}\) 


        

<hr>

------------------------------


