__author__ = 'veronika'
from django import template

register = template.Library()


@register.filter(name='find_val')
def find_val(value, arg):
    return value[arg]


@register.filter(name='remove_def')
def remove_def(value):
    return value.keys()


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
    for i in clist:
        color_p += i
    return color_p
