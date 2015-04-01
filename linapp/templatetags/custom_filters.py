__author__ = 'veronika'
from django import template

register = template.Library()


@register.filter(name='find_val')
def find_val(value, arg):
    return value[arg]


@register.filter(name='get_keys')
def get_keys(value):
    return sorted(value.keys())


@register.filter(name='has_value')
def has_value(value, arg):
    if value[arg]:
        return True
    else:
        return False


@register.filter(name='color')
def color(value, arg):
    color_p = ''
    clist = value[arg]
    print clist
    for i in clist:
        color_p += i
    return color_p
